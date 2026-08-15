"""Independent verifier for structured translated-copy CEGIS artifacts.

This script does not import the search implementation.  It reconstructs every
listed edge from the base graph, the author-supplied Heule template, and the
two exact anchors of each translated copy.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from exact_geometry import Point, parse_mathematica_vertices
from graph_verifier import load_graph_json


def _parse_edge_file(path: str | Path) -> tuple[int, set[tuple[int, int]]]:
    vertex_count = declared_edges = None
    edges: set[tuple[int, int]] = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.split()
            if not fields or fields[0] == "c":
                continue
            if fields[0] == "p":
                if len(fields) != 4 or fields[1] != "edge":
                    raise ValueError(f"{path}:{line_number}: invalid header")
                vertex_count, declared_edges = int(fields[2]), int(fields[3])
            elif fields[0] == "e" and len(fields) == 3:
                left, right = int(fields[1]) - 1, int(fields[2]) - 1
                edges.add((min(left, right), max(left, right)))
            else:
                raise ValueError(f"{path}:{line_number}: invalid record")
    if vertex_count is None or declared_edges != len(edges):
        raise ValueError(f"{path}: edge header mismatch")
    return vertex_count, edges


def _circle_translation(r: Point) -> tuple[Point, sp.Poly, bool]:
    d2 = sp.radsimp(sp.expand(r.x * r.x + r.y * r.y))
    h_squared = sp.radsimp(sp.cancel((4 - d2) / (4 * d2)))
    h = sp.sqrt(h_squared)
    translation = Point(
        r.x / 2 + r.y * h,
        r.y / 2 - r.x * h,
    )
    variable = sp.symbols("x")
    polynomial = sp.Poly(
        sp.minimal_polynomial(h, variable),
        variable,
        domain=sp.QQ,
    )
    group, _ = polynomial.galois_group()
    return translation, polynomial, not bool(group.is_abelian)


def verify_translated_copy_artifact(
    graph_path: str | Path,
    report_path: str | Path,
    base_graph_path: str | Path,
    source_vertex_path: str | Path,
    source_edge_path: str | Path,
) -> dict[str, object]:
    """Reconstruct coordinates, listed edges, and the saved five-coloring."""

    points, listed_edges = load_graph_json(graph_path)
    base_points, base_edges = load_graph_json(base_graph_path)
    source_points = parse_mathematica_vertices(source_vertex_path)
    source_count, source_edges = _parse_edge_file(source_edge_path)
    with Path(report_path).open(encoding="utf-8") as handle:
        report = json.load(handle)
    selected = [
        int(index)
        for index in report["selected_source_r_indices_zero_based"]
    ]
    coloring_raw = report.get("whole_graph_coloring")
    coloring = (
        tuple(int(color) for color in coloring_raw)
        if coloring_raw is not None
        else None
    )

    if source_count != len(source_points):
        raise ValueError("source coordinate and edge counts disagree")
    expected_vertex_count = len(base_points) + len(selected) * source_count
    coordinates_match = len(points) == expected_vertex_count
    if coordinates_match:
        coordinates_match = points[: len(base_points)] == base_points

    # Verify each base and template edge once.  Translation cancels from every
    # copied template edge, so this symbolically verifies all their instances.
    base_edges_are_unit = all(
        base_points[left].squared_distance(base_points[right]) == 1
        for left, right in base_edges
    )
    source_edges_are_unit = all(
        source_points[left].squared_distance(source_points[right]) == 1
        for left, right in source_edges
    )

    lookup = {point: index for index, point in enumerate(base_points)}
    expected_edges = set(base_edges)
    translation_certificates: list[dict[str, object]] = []
    if len(lookup) != len(base_points):
        raise ValueError("base graph has duplicate exact coordinates")

    for copy_index, r_index in enumerate(selected):
        r = source_points[r_index]
        translation, polynomial, nonabelian = _circle_translation(r)
        origin_unit = translation.squared_distance(Point(0, 0)) == 1
        second_anchor_unit = translation.squared_distance(r) == 1
        offset = len(base_points) + copy_index * source_count
        expected_block = [
            Point(point.x + translation.x, point.y + translation.y)
            for point in source_points
        ]
        if coordinates_match:
            coordinates_match = (
                points[offset : offset + source_count] == expected_block
            )
        for left, right in source_edges:
            expected_edges.add((offset + left, offset + right))
        cross_count = 0
        for source_index, source_point in enumerate(source_points):
            anchors = (
                source_point,
                Point(source_point.x + r.x, source_point.y + r.y),
            )
            for anchor in anchors:
                base_index = lookup.get(anchor)
                if base_index is not None:
                    expected_edges.add((base_index, offset + source_index))
                    cross_count += 1
        translation_certificates.append(
            {
                "source_r_index_zero_based": r_index,
                "minimal_polynomial": sp.sstr(polynomial.as_expr()),
                "galois_group_nonabelian": nonabelian,
                "origin_anchor_squared_distance_one": origin_unit,
                "second_anchor_squared_distance_one": second_anchor_unit,
                "cross_edges": cross_count,
            }
        )

    coloring_conflicts = (
        None
        if coloring is None or len(coloring) != len(points)
        else sum(coloring[left] == coloring[right] for left, right in listed_edges)
    )
    duplicate_vertices = len(points) - len(set(points))
    edges_match = expected_edges == listed_edges
    translation_certificates_valid = all(
        certificate["galois_group_nonabelian"]
        and certificate["origin_anchor_squared_distance_one"]
        and certificate["second_anchor_squared_distance_one"]
        for certificate in translation_certificates
    )
    valid = bool(
        coordinates_match
        and base_edges_are_unit
        and source_edges_are_unit
        and translation_certificates_valid
        and duplicate_vertices == 0
        and edges_match
        and coloring_conflicts == 0
    )
    return {
        "valid": valid,
        "vertices": len(points),
        "listed_edges": len(listed_edges),
        "duplicate_vertices": duplicate_vertices,
        "coordinates_match_reconstruction": coordinates_match,
        "base_edges_symbolically_unit": base_edges_are_unit,
        "source_template_edges_symbolically_unit": source_edges_are_unit,
        "edges_match_reconstruction": edges_match,
        "missing_listed_edges": len(expected_edges - listed_edges),
        "unexpected_listed_edges": len(listed_edges - expected_edges),
        "coloring_conflicts": coloring_conflicts,
        "translations": translation_certificates,
        "scope_note": (
            "This verifies the exact listed-edge subgraph.  The search artifact "
            "intentionally does not claim completeness for unit edges between "
            "different translated copies."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--base-graph", type=Path, required=True)
    parser.add_argument("--source-vtx", type=Path, required=True)
    parser.add_argument("--source-edge", type=Path, required=True)
    args = parser.parse_args()
    result = verify_translated_copy_artifact(
        args.graph,
        args.report,
        args.base_graph,
        args.source_vtx,
        args.source_edge,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
