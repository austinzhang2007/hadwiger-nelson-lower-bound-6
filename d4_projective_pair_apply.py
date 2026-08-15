"""Apply an exact all-projective forced-pair scan to a D4 CEGIS graph."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Mapping, Sequence

import sympy as sp

from coloring_sat import ColoringSAT, parse_dimacs_edge
from d4_blocker_cegis import (
    _projective_from_payload,
    _projective_payload,
    _write_current_artifacts,
    _write_json,
)
from d4_exact import D4ProjectivePoint, primary_d4_coordinate_backend
from exact_geometry import parse_mathematica_vertices
from graph_verifier import load_graph_json


def validate_pair_scan(
    payload: Mapping[str, object],
) -> tuple[dict[int, Mapping[str, object]], tuple[tuple[int, int], ...]]:
    """Validate sparse scan references before touching a graph."""

    raw_records = payload.get("used_candidates")
    raw_edges = payload.get("pair_edges")
    if not isinstance(raw_records, list) or not isinstance(raw_edges, list):
        raise ValueError("pair scan needs candidate and edge lists")
    records: dict[int, Mapping[str, object]] = {}
    for raw in raw_records:
        if not isinstance(raw, Mapping) or not isinstance(
            raw.get("candidate_index"), int
        ):
            raise ValueError("invalid candidate record")
        index = int(raw["candidate_index"])
        if index in records:
            raise ValueError("duplicate candidate index")
        records[index] = raw
    edges: set[tuple[int, int]] = set()
    for raw in raw_edges:
        if (
            not isinstance(raw, list)
            or len(raw) != 2
            or not all(isinstance(value, int) for value in raw)
        ):
            raise ValueError("invalid pair edge")
        left, right = sorted((int(raw[0]), int(raw[1])))
        if left == right:
            raise ValueError("pair edge is a loop")
        if left not in records or right not in records:
            raise ValueError("pair edge references a missing candidate")
        edges.add((left, right))
    return records, tuple(sorted(edges))


def apply_projective_forced_pair_scan(
    *,
    base_graph_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    resume_report_path: str | Path,
    scan_path: str | Path,
    output_directory: str | Path,
    symmetry_clique: Sequence[int] = (548, 1149, 668),
) -> dict[str, object]:
    """Add exact forced candidates and their certified internal unit edges."""

    started = time.monotonic()
    base_points, base_edges = load_graph_json(base_graph_path)
    with Path(primary_report_path).open(encoding="utf-8") as handle:
        primary = json.load(handle)
    selected = [
        int(index) for index in primary["selected_source_r_indices_zero_based"]
    ]
    source_points = parse_mathematica_vertices(source_vertex_path)
    original_base_count = len(base_points) - len(selected) * len(source_points)
    affine = primary_d4_coordinate_backend(
        base_points[:original_base_count], source_points, selected
    )
    with Path(resume_report_path).open(encoding="utf-8") as handle:
        resume = json.load(handle)
    with Path(scan_path).open(encoding="utf-8") as handle:
        scan = json.load(handle)
    records_by_index, pair_edges = validate_pair_scan(scan)

    vertex_count, current_edges = parse_dimacs_edge(resume["edge_path"])
    iterations = list(resume["iterations"])
    candidate_projective: list[D4ProjectivePoint] = [
        affine.algebra.normalize_projective(
            _projective_from_payload(
                {
                    **record["projective_certificate"],
                    "fifth_neighbor_zero_based": record[
                        "projective_certificate"
                    ].get(
                        "fifth_neighbor_zero_based",
                        record.get("generator_4_zero_based", -1),
                    ),
                }
            )
        )
        for record in iterations
    ]
    candidate_points = [affine.materialize(point) for point in candidate_projective]
    if vertex_count != len(base_points) + len(candidate_points):
        raise ValueError("resume vertex and certificate counts differ")
    coloring = tuple(int(color) for color in resume["whole_graph_coloring"])
    if not ColoringSAT(
        vertex_count, current_edges, 5, break_color_symmetry=False
    ).validate(coloring):
        raise ValueError("resume coloring is invalid")

    exact = affine.as_projective_backend(candidate_projective)
    all_points = [*base_points, *candidate_points]
    by_key: dict[tuple[float, float], list[int]] = {}
    for vertex, point in enumerate(all_points):
        x, y = point.approximate()
        by_key.setdefault((round(x, 10), round(y, 10)), []).append(vertex)

    def projective_at(vertex: int) -> D4ProjectivePoint:
        if vertex < len(exact.coordinates):
            return exact.coordinates[vertex]
        return candidate_projective[vertex - len(base_points)]

    scan_to_vertex: dict[int, int] = {}
    added_vertices = 0
    added_generator_edges = 0
    missing_color = int(scan["missing_color"])
    present_colors = set(range(5)) - {missing_color}
    for scan_index, record in sorted(records_by_index.items()):
        certificate = record.get("projective_certificate")
        neighbors_raw = record.get("neighbors_zero_based")
        if not isinstance(certificate, Mapping) or not isinstance(
            neighbors_raw, list
        ):
            raise ValueError("candidate lacks exact certificate or neighbors")
        neighbors = tuple(sorted({int(index) for index in neighbors_raw}))
        if any(index < 0 or index >= vertex_count for index in neighbors):
            raise ValueError("candidate neighbor outside resume graph")
        if {coloring[index] for index in neighbors} != present_colors:
            raise ValueError("candidate is not forced to the declared color")
        projective = affine.algebra.normalize_projective(
            _projective_from_payload(dict(certificate))
        )
        point = affine.materialize(projective)
        x, y = point.approximate()
        key = (round(x, 10), round(y, 10))
        vertex = next(
            (
                known
                for known in by_key.get(key, ())
                if exact.projective_equal(projective_at(known), projective)
            ),
            -1,
        )
        if vertex < 0:
            vertex = len(base_points) + len(candidate_points)
            candidate_points.append(point)
            candidate_projective.append(projective)
            by_key.setdefault(key, []).append(vertex)
            iterations.append(
                {
                    "iteration": len(iterations) + 1,
                    "projective_forced_pair_scan": str(scan_path),
                    "candidate_vertex_zero_based": vertex,
                    "x_exact": sp.sstr(point.x),
                    "y_exact": sp.sstr(point.y),
                    "forced_color": missing_color,
                    "listed_neighbors_zero_based": list(neighbors),
                    "projective_certificate": {
                        **_projective_payload(projective),
                        "fifth_neighbor_zero_based": -1,
                    },
                }
            )
            added_vertices += 1
        scan_to_vertex[scan_index] = vertex
        for neighbor in neighbors:
            edge = (min(neighbor, vertex), max(neighbor, vertex))
            if edge not in current_edges:
                current_edges.add(edge)
                added_generator_edges += 1

    added_pair_edges = 0
    for left, right in pair_edges:
        left_vertex = scan_to_vertex[left]
        right_vertex = scan_to_vertex[right]
        if not exact.projective_unit(
            projective_at(left_vertex), projective_at(right_vertex)
        ):
            raise ValueError("scan contains a non-unit forced pair")
        edge = (min(left_vertex, right_vertex), max(left_vertex, right_vertex))
        if edge not in current_edges:
            current_edges.add(edge)
            added_pair_edges += 1

    output_directory = Path(output_directory)
    checkpoint = _write_current_artifacts(
        output_directory=output_directory,
        base_points=base_points,
        base_edges=base_edges,
        candidate_points=candidate_points,
        current_edges=current_edges,
        iterations=iterations,
        coloring=None,
        symmetry_clique=symmetry_clique,
        status="PROJECTIVE_FORCED_PAIRS_PENDING_SAT",
        complete_unit_distance_graph=False,
    )
    phases = [
        vertex * 5 + color + 1
        if color == coloring[vertex]
        else -(vertex * 5 + color + 1)
        for vertex in range(len(coloring))
        for color in range(5)
    ]
    solve_started = time.monotonic()
    next_coloring = ColoringSAT(
        len(base_points) + len(candidate_points),
        current_edges,
        5,
        break_color_symmetry=False,
        symmetry_clique=symmetry_clique,
    ).solve(phases=phases)
    solve_seconds = time.monotonic() - solve_started
    status = (
        "PROJECTIVE_FORCED_PAIR_GRAPH_5_UNSAT_WITHOUT_PROOF"
        if next_coloring is None
        else "PROJECTIVE_FORCED_PAIRS_WITH_5_COLORING"
    )
    payload = _write_current_artifacts(
        output_directory=output_directory,
        base_points=base_points,
        base_edges=base_edges,
        candidate_points=candidate_points,
        current_edges=current_edges,
        iterations=iterations,
        coloring=next_coloring,
        symmetry_clique=symmetry_clique,
        status=status,
        complete_unit_distance_graph=False,
    )
    payload["projective_forced_pair_application"] = {
        "source_report": str(resume_report_path),
        "scan_path": str(scan_path),
        "missing_color": missing_color,
        "scan_pair_edges": len(pair_edges),
        "added_vertices": added_vertices,
        "added_generator_edges": added_generator_edges,
        "added_pair_edges": added_pair_edges,
        "solve_seconds": solve_seconds,
    }
    payload["elapsed_seconds"] = time.monotonic() - started
    report_path = Path(str(payload["graph_path"])).with_suffix("").with_suffix(
        ".json"
    )
    _write_json(report_path, payload)
    return payload
