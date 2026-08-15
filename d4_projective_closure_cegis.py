"""CEGIS that closes five-color blockers over previously added D4 points.

Unlike :mod:`d4_blocker_cegis`, every already certified candidate may be used
as a unit-circle center.  Coordinates are kept projective in the fixed
32-dimensional D4 algebra, so all retained incidences are exact.
"""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import time
from typing import Callable, Sequence

import sympy as sp

from blocking_point_search import exact_five_color_blockers
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


def _report_path(payload: dict[str, object]) -> Path:
    return Path(str(payload["graph_path"])).with_suffix("").with_suffix(".json")


def run_d4_projective_closure_cegis(
    *,
    base_graph_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    resume_report_path: str | Path,
    output_directory: str | Path,
    max_rounds: int = 5,
    blockers_per_round: int = 25,
    symmetry_clique: Sequence[int] = (548, 1149, 668),
    progress: Callable[[str, dict[str, object]], None] | None = None,
) -> dict[str, object]:
    """Find exact blockers using all current vertices as circle centers."""

    if max_rounds < 1 or blockers_per_round < 1:
        raise ValueError("round and blocker limits must be positive")
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
    affine_backend = primary_d4_coordinate_backend(
        base_points[:original_base_count],
        source_points,
        selected,
    )

    with Path(resume_report_path).open(encoding="utf-8") as handle:
        resume = json.load(handle)
    vertex_count, current_edges = parse_dimacs_edge(resume["edge_path"])
    iterations = list(resume["iterations"])
    candidate_projective: list[D4ProjectivePoint] = [
        affine_backend.algebra.normalize_projective(
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
    candidate_points = [
        affine_backend.materialize(projective)
        for projective in candidate_projective
    ]
    if vertex_count != len(base_points) + len(candidate_points):
        raise ValueError("resume vertex and certificate counts differ")
    coloring = tuple(int(color) for color in resume["whole_graph_coloring"])
    if not ColoringSAT(
        vertex_count,
        current_edges,
        5,
        break_color_symmetry=False,
    ).validate(coloring):
        raise ValueError("resume coloring is invalid")

    by_key: dict[tuple[float, float], list[int]] = {}
    for index, point in enumerate(candidate_points):
        x, y = point.approximate()
        by_key.setdefault((round(x, 10), round(y, 10)), []).append(index)
    closure_rounds = list(resume.get("projective_closure_rounds", []))
    latest: dict[str, object] = {}
    status = "PROJECTIVE_CLOSURE_ITERATION_LIMIT_WITH_5_COLORING"

    for local_round in range(1, max_rounds + 1):
        all_points = [*base_points, *candidate_points]
        exact_backend = affine_backend.as_projective_backend(candidate_projective)
        scan_started = time.monotonic()
        blockers, stats = exact_five_color_blockers(
            all_points,
            coloring,
            exact_backend=exact_backend,
            max_blockers=blockers_per_round,
            progress=(
                (lambda event, data: progress(event, {"round": local_round, **data}))
                if progress is not None
                else None
            ),
        )
        scan_seconds = time.monotonic() - scan_started
        if not blockers:
            status = "PROJECTIVE_SINGLE_POINT_POOL_STAGNATED_WITH_5_COLORING"
            break

        edges_before = len(current_edges)
        vertices_before = len(candidate_points)
        records: list[dict[str, object]] = []
        for blocker in blockers:
            neighbor_by_color = {
                int(coloring[index]): index for index in blocker.neighbors
            }
            projective = exact_backend.recover_projective(
                blocker.generator_pair,
                (neighbor_by_color[2], neighbor_by_color[3]),
                (neighbor_by_color[4],),
            )
            if projective is None:
                raise AssertionError("closure blocker lost its exact certificate")
            x, y = blocker.point.approximate()
            key = (round(x, 10), round(y, 10))
            if any(
                affine_backend.projective_equal(
                    candidate_projective[index],
                    projective,
                )
                for index in by_key.get(key, ())
            ):
                continue
            vertex = len(base_points) + len(candidate_points)
            current_edges.update(
                (min(neighbor, vertex), max(neighbor, vertex))
                for neighbor in blocker.neighbors
            )
            candidate_points.append(blocker.point)
            candidate_projective.append(projective)
            by_key.setdefault(key, []).append(len(candidate_points) - 1)
            records.append(
                {
                    "iteration": len(iterations) + len(records) + 1,
                    "projective_closure_round": len(closure_rounds) + 1,
                    "candidate_vertex_zero_based": vertex,
                    "x_exact": sp.sstr(blocker.point.x),
                    "y_exact": sp.sstr(blocker.point.y),
                    "generator_01_zero_based": list(blocker.generator_pair),
                    "generator_23_zero_based": [
                        neighbor_by_color[2],
                        neighbor_by_color[3],
                    ],
                    "generator_4_zero_based": neighbor_by_color[4],
                    "listed_neighbors_zero_based": list(blocker.neighbors),
                    "projective_certificate": {
                        **_projective_payload(projective),
                        "fifth_neighbor_zero_based": projective.fifth_neighbor,
                    },
                }
            )
        if not records:
            status = "PROJECTIVE_SINGLE_POINT_POOL_REPEATED_WITH_5_COLORING"
            break
        iterations.extend(records)
        if progress is not None:
            progress(
                "closure_blockers_added",
                {
                    "round": local_round,
                    "new_vertices": len(candidate_points) - vertices_before,
                    "new_edges": len(current_edges) - edges_before,
                    "scan_seconds": scan_seconds,
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
        closure_rounds.append(
            {
                "round": len(closure_rounds) + 1,
                "blocking_search": asdict(stats),
                "new_vertices": len(candidate_points) - vertices_before,
                "new_edges": len(current_edges) - edges_before,
                "scan_seconds": scan_seconds,
                "solve_seconds": solve_seconds,
                "five_sat_after_add": next_coloring is not None,
            }
        )
        if progress is not None:
            progress(
                "sat_result",
                {
                    "round": local_round,
                    "five_sat": next_coloring is not None,
                    "solve_seconds": solve_seconds,
                    "candidate_count": len(candidate_points),
                },
            )
        if next_coloring is None:
            coloring = ()
            status = "EXACT_LISTED_PROJECTIVE_CLOSURE_5_UNSAT_WITHOUT_PROOF"
        else:
            coloring = next_coloring
            status = "PROJECTIVE_CLOSURE_RUNNING_WITH_5_COLORING"
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
        latest["projective_closure_rounds"] = closure_rounds
        latest["elapsed_seconds"] = time.monotonic() - started
        _write_json(_report_path(latest), latest)
        if next_coloring is None:
            break

    if status == "PROJECTIVE_CLOSURE_RUNNING_WITH_5_COLORING":
        status = "PROJECTIVE_CLOSURE_ITERATION_LIMIT_WITH_5_COLORING"
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
        latest["projective_closure_rounds"] = closure_rounds
        latest["elapsed_seconds"] = time.monotonic() - started
        _write_json(_report_path(latest), latest)
    return latest
