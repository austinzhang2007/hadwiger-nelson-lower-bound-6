"""CEGIS using adjacent pairs of four-color-forced D4 candidates."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Callable, Sequence

import sympy as sp

from coloring_sat import ColoringSAT, parse_dimacs_edge
from d4_blocker_cegis import (
    _projective_from_payload,
    _projective_payload,
    _write_current_artifacts,
    _write_json,
)
from d4_exact import D4ProjectivePoint, primary_d4_coordinate_backend
from d4_pair_search import exact_forced_pair_edges, exact_four_color_forced_points
from exact_geometry import Point, parse_mathematica_vertices
from graph_verifier import load_graph_json


def _report_path(payload: dict[str, object]) -> Path:
    return Path(str(payload["graph_path"])).with_suffix("").with_suffix(".json")


def run_d4_pair_cegis(
    *,
    base_graph_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    resume_report_path: str | Path,
    output_directory: str | Path,
    max_rounds: int = 5,
    max_pairs_per_round: int | None = None,
    symmetry_clique: Sequence[int] = (548, 1149, 668),
    progress: Callable[[str, dict[str, object]], None] | None = None,
) -> dict[str, object]:
    """Add every selected exact forced pair, then solve the listed-edge graph."""

    started = time.monotonic()
    base_points, base_edges = load_graph_json(base_graph_path)
    with Path(primary_report_path).open(encoding="utf-8") as handle:
        primary_report = json.load(handle)
    selected = [
        int(index)
        for index in primary_report["selected_source_r_indices_zero_based"]
    ]
    source_points = parse_mathematica_vertices(source_vertex_path)
    original_base_count = len(base_points) - len(selected) * len(source_points)
    backend = primary_d4_coordinate_backend(
        base_points[:original_base_count],
        source_points,
        selected,
    )

    with Path(resume_report_path).open(encoding="utf-8") as handle:
        resume = json.load(handle)
    vertex_count, current_edges = parse_dimacs_edge(resume["edge_path"])
    iterations = list(resume["iterations"])
    candidate_projective: list[D4ProjectivePoint] = []
    for record in iterations:
        certificate = record["projective_certificate"]
        candidate_projective.append(
            backend.algebra.normalize_projective(_projective_from_payload(
                {
                    **certificate,
                    "fifth_neighbor_zero_based": certificate.get(
                        "fifth_neighbor_zero_based",
                        record.get("generator_4_zero_based", -1),
                    ),
                }
            ))
        )
    candidate_points = [
        backend.materialize(projective)
        for projective in candidate_projective
    ]
    if vertex_count != len(base_points) + len(candidate_points):
        raise ValueError("resume vertex and candidate counts differ")
    coloring = tuple(int(color) for color in resume["whole_graph_coloring"])
    if not ColoringSAT(
        vertex_count,
        current_edges,
        5,
        break_color_symmetry=False,
    ).validate(coloring):
        raise ValueError("resume coloring is invalid")
    pair_rounds = list(resume.get("pair_cegis_rounds", []))

    by_key: dict[tuple[float, float], list[int]] = {}
    for index, point in enumerate(candidate_points):
        x, y = point.approximate()
        by_key.setdefault((round(x, 10), round(y, 10)), []).append(index)

    latest: dict[str, object] = {}
    status = "PAIR_ITERATION_LIMIT_WITH_5_COLORING"
    for round_number in range(1, max_rounds + 1):
        base_coloring = coloring[: len(base_points)]
        forced_groups = []
        pair_groups = []
        search_counts: list[dict[str, int]] = []
        for missing_color in range(5):
            forced, stats = exact_four_color_forced_points(
                base_points,
                base_coloring,
                backend,
                missing_color=missing_color,
            )
            pair_edges = sorted(exact_forced_pair_edges(forced, backend))
            forced_groups.append(forced)
            pair_groups.append(pair_edges)
            search_counts.append(
                {
                    "missing_color": missing_color,
                    "forced_candidates": len(forced),
                    "unit_pairs": len(pair_edges),
                    "common_intersections": stats.common_intersections,
                }
            )
        available_pairs = [
            (missing, left, right)
            for missing, edges in enumerate(pair_groups)
            for left, right in edges
        ]
        if max_pairs_per_round is not None:
            available_pairs = available_pairs[:max_pairs_per_round]
        if not available_pairs:
            status = "FORCED_PAIR_POOL_STAGNATED_WITH_5_COLORING"
            break

        edges_before = len(current_edges)
        vertices_before = len(candidate_points)
        pair_edges_added = 0

        def ensure_candidate(candidate) -> int:
            x, y = candidate.point.approximate()
            key = (round(x, 10), round(y, 10))
            for candidate_index in by_key.get(key, ()):
                if backend.projective_equal(
                    candidate_projective[candidate_index],
                    candidate.projective,
                ):
                    vertex = len(base_points) + candidate_index
                    current_edges.update(
                        (min(neighbor, vertex), max(neighbor, vertex))
                        for neighbor in candidate.neighbors
                    )
                    return vertex

            candidate_index = len(candidate_points)
            vertex = len(base_points) + candidate_index
            candidate_points.append(candidate.point)
            candidate_projective.append(candidate.projective)
            by_key.setdefault(key, []).append(candidate_index)
            current_edges.update(
                (min(neighbor, vertex), max(neighbor, vertex))
                for neighbor in candidate.neighbors
            )
            iterations.append(
                {
                    "iteration": len(iterations) + 1,
                    "pair_cegis_round": round_number,
                    "candidate_vertex_zero_based": vertex,
                    "x_exact": sp.sstr(candidate.point.x),
                    "y_exact": sp.sstr(candidate.point.y),
                    "forced_color": candidate.forced_color,
                    "listed_base_neighbors_zero_based": list(candidate.neighbors),
                    "projective_certificate": {
                        **_projective_payload(candidate.projective),
                        "fifth_neighbor_zero_based": -1,
                    },
                }
            )
            return vertex

        for missing, left, right in available_pairs:
            candidates = forced_groups[missing]
            left_vertex = ensure_candidate(candidates[left])
            right_vertex = ensure_candidate(candidates[right])
            edge = (min(left_vertex, right_vertex), max(left_vertex, right_vertex))
            if edge not in current_edges:
                current_edges.add(edge)
                pair_edges_added += 1

        if len(current_edges) == edges_before:
            status = "FORCED_PAIR_POOL_REPEATED_WITH_5_COLORING"
            break
        if progress is not None:
            progress(
                "pair_batch_added",
                {
                    "round": round_number,
                    "available_pairs": len(available_pairs),
                    "new_vertices": len(candidate_points) - vertices_before,
                    "new_edges": len(current_edges) - edges_before,
                    "new_pair_edges": pair_edges_added,
                },
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
        pair_rounds.append(
            {
                "round": len(pair_rounds) + 1,
                "search_counts": search_counts,
                "selected_pairs": len(available_pairs),
                "new_vertices": len(candidate_points) - vertices_before,
                "new_edges": len(current_edges) - edges_before,
                "new_pair_edges": pair_edges_added,
                "five_sat_after_add": next_coloring is not None,
                "solve_seconds": solve_seconds,
            }
        )
        if progress is not None:
            progress(
                "sat_result",
                {
                    "round": round_number,
                    "five_sat": next_coloring is not None,
                    "solve_seconds": solve_seconds,
                    "candidate_count": len(candidate_points),
                },
            )
        if next_coloring is None:
            coloring = ()
            status = "EXACT_LISTED_PAIR_GRAPH_5_UNSAT_WITHOUT_PROOF"
        else:
            coloring = next_coloring
            status = "PAIR_RUNNING_WITH_5_COLORING"
        latest = _write_current_artifacts(
            output_directory=Path(output_directory),
            base_points=base_points,
            base_edges=base_edges,
            candidate_points=candidate_points,
            current_edges=current_edges,
            iterations=iterations,
            coloring=coloring if coloring else None,
            symmetry_clique=symmetry_clique,
            status=status,
            complete_unit_distance_graph=False,
        )
        latest["pair_cegis_rounds"] = pair_rounds
        latest["elapsed_seconds"] = time.monotonic() - started
        _write_json(_report_path(latest), latest)
        if next_coloring is None:
            break

    if status == "PAIR_RUNNING_WITH_5_COLORING":
        status = "PAIR_ITERATION_LIMIT_WITH_5_COLORING"
    if not latest or latest.get("status") != status:
        latest = _write_current_artifacts(
            output_directory=Path(output_directory),
            base_points=base_points,
            base_edges=base_edges,
            candidate_points=candidate_points,
            current_edges=current_edges,
            iterations=iterations,
            coloring=coloring,
            symmetry_clique=symmetry_clique,
            status=status,
            complete_unit_distance_graph=False,
        )
        latest["pair_cegis_rounds"] = pair_rounds
        latest["elapsed_seconds"] = time.monotonic() - started
        _write_json(_report_path(latest), latest)
    return latest
