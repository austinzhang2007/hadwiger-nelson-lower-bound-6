"""Persisted exact listed-edge audit with independent GMP arithmetic.

Coordinates come from the original base prefix, author translation recipe,
and homogeneous certificates. Redundant radical JSON consistency, duplicate
absence, edge completeness and UNSAT proofs remain separate gates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import ijson
from primary_field_audit import PrimaryField


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def audit_coordinates(coordinates, edges, coloring, *, output,
                      batch_size=1000, progress=None, metadata=None):
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    n = len(coordinates)
    f = PrimaryField()
    if len(coloring) != n or any(type(c) is not int or not 0 <= c < 5
                                 for c in coloring):
        raise ValueError("invalid five-color array")
    if any(len(p) != 3 or any(len(v) != 32 for v in p) for p in coordinates):
        raise ValueError("invalid coordinate dimensions")
    if any(p[2] == f.zero for p in coordinates):
        raise ValueError("zero projective denominator")
    ordered = sorted(set(tuple(edge) for edge in edges))
    if any(not (0 <= u < v < n) for u, v in ordered):
        raise ValueError("invalid edge")
    started = time.monotonic()
    invalid = []
    conflicts = sum(coloring[u] == coloring[v] for u, v in ordered)
    common = {
        "metadata": metadata or {}, "vertices": n, "edges": len(ordered),
        "class_sizes": [coloring.count(c) for c in range(5)],
        "coloring_conflicts": conflicts, "coloring_valid": conflicts == 0,
        "coloring_sha256": hashlib.sha256(json.dumps(
            coloring, separators=(",", ":")).encode()).hexdigest(),
        "complete_unit_graph_certified": False,
        "duplicate_free_certified": False,
        "redundant_graph_json_consistency_certified": False,
        "lower_bound_6_proved": False,
    }

    def save(checked, finished=False):
        payload = {
            **common, "status": "FINISHED" if finished else "RUNNING",
            "checked_edges": checked, "invalid_edges": invalid,
            "listed_edges_exactly_unit": (not invalid) if finished else None,
            "elapsed_seconds": time.monotonic() - started,
        }
        write_json(output, payload)
        if progress:
            progress(payload)
        return payload

    for checked, (u, v) in enumerate(ordered, 1):
        if not f.unit(coordinates[u], coordinates[v]):
            invalid.append([u, v])
        if checked % batch_size == 0:
            save(checked)
    return save(len(ordered), finished=True)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base-graph", "primary-report", "source-vertices", "report", "edge", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--batch-size", type=int, default=1000)
    args = parser.parse_args()
    started = time.monotonic()

    def event(stage, **details):
        payload = {"status": "LOADING", "stage": stage,
                   "elapsed_seconds": time.monotonic() - started, **details}
        write_json(args.output, payload)
        print(json.dumps(payload, sort_keys=True), flush=True)

    event("author_base")
    # Reconstruction only; the edge arithmetic imports no search code.
    from d4_exact import primary_d4_coordinate_backend
    from exact_geometry import Point, parse_mathematica_vertices, parse_real_expr
    from graph_verifier import _parse_dimacs_edge

    primary = json.loads(args.primary_report.read_text())
    selected = primary["selected_source_r_indices_zero_based"]
    source = parse_mathematica_vertices(args.source_vertices)
    base_count = int(primary["vertices"])
    original_count = base_count - len(selected) * len(source)
    if original_count <= 0:
        raise ValueError("invalid author-copy recipe")
    original = []
    with args.base_graph.open("rb") as stream:
        for index, vertex in enumerate(ijson.items(stream, "vertices.item")):
            if index == original_count:
                break
            if vertex["id"] != index + 1:
                raise ValueError("nonconsecutive base vertex id")
            original.append(Point(parse_real_expr(vertex["x"]), parse_real_expr(vertex["y"])))
    if len(original) != original_count:
        raise ValueError("base prefix shorter than recipe")
    affine = primary_d4_coordinate_backend(original, source, selected)
    field = PrimaryField()
    coordinates = [(field.element(x), field.element(y), field.one)
                   for x, y in affine.coordinates]
    if len(coordinates) != base_count:
        raise ValueError("base recipe count mismatch")
    event("coefficient_certificates", base_vertices=base_count)
    with args.report.open("rb") as stream:
        for record in ijson.items(stream, "iterations.item"):
            cert = record["projective_certificate"]
            if cert.get("basis_masks") != list(range(32)):
                raise ValueError("unsupported field basis")
            coordinates.append(tuple(field.element(cert[key]) for key in (
                "numerator_x_coefficients", "numerator_y_coefficients", "denominator_coefficients")))
    with args.report.open("rb") as stream:
        coloring = list(ijson.items(stream, "whole_graph_coloring.item"))
    n, edges = _parse_dimacs_edge(args.edge)
    if n != len(coordinates):
        raise ValueError("coordinate and edge vertex counts differ")
    event("input_hashes", vertices=n, edges=len(edges))
    metadata = {
        "coordinate_definition": "original base prefix + author translations + homogeneous certificates",
        "input_sha256": {str(path): sha256(path) for path in (
            args.base_graph, args.primary_report, args.source_vertices, args.report, args.edge)},
        "arithmetic": "primary_field_audit.PrimaryField / GMP exact rationals",
        "arithmetic_sha256": sha256(Path(__file__).with_name("primary_field_audit.py")),
    }

    def progress(payload):
        print(json.dumps({key: payload[key] for key in (
            "status", "checked_edges", "edges", "elapsed_seconds",
            "listed_edges_exactly_unit", "coloring_conflicts")}), flush=True)

    result = audit_coordinates(coordinates, edges, coloring, output=args.output,
                               batch_size=args.batch_size, progress=progress, metadata=metadata)
    return 0 if result["listed_edges_exactly_unit"] and result["coloring_valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
