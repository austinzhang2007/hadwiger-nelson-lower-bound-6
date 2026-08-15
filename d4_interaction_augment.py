"""Add every numerically located and exactly certified D4 interaction edge."""

from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Callable, Sequence

import numpy as np
from scipy.spatial import cKDTree

from coloring_sat import ColoringSAT, parse_dimacs_edge
from d4_blocker_cegis import (
    _projective_from_payload,
    _write_current_artifacts,
    _write_json,
)
from d4_exact import primary_d4_coordinate_backend
from exact_geometry import parse_mathematica_vertices
from graph_verifier import load_graph_json


def augment_exact_interactions(
    *,
    base_graph_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    resume_report_path: str | Path,
    output_directory: str | Path,
    tolerance: float = 1e-8,
    symmetry_clique: Sequence[int] = (548, 1149, 668),
    progress: Callable[[str, dict[str, object]], None] | None = None,
) -> dict[str, object]:
    """Numerically locate unit pairs, certify them exactly, and re-solve."""

    started = time.monotonic()
    base_points, base_edges = load_graph_json(base_graph_path)
    with Path(primary_report_path).open(encoding="utf-8") as handle:
        primary = json.load(handle)
    selected = [
        int(index)
        for index in primary["selected_source_r_indices_zero_based"]
    ]
    source_points = parse_mathematica_vertices(source_vertex_path)
    original_base_count = len(base_points) - len(selected) * len(source_points)
    affine = primary_d4_coordinate_backend(
        base_points[:original_base_count],
        source_points,
        selected,
    )
    with Path(resume_report_path).open(encoding="utf-8") as handle:
        resume = json.load(handle)
    iterations = list(resume["iterations"])
    candidates = [
        affine.algebra.normalize_projective(
            _projective_from_payload(
                {
                    **record["projective_certificate"],
                    "fifth_neighbor_zero_based": record[
                        "projective_certificate"
                    ].get("fifth_neighbor_zero_based", -1),
                }
            )
        )
        for record in iterations
    ]
    candidate_points = [affine.materialize(point) for point in candidates]
    points = [*base_points, *candidate_points]
    exact = affine.as_projective_backend(candidates)
    vertex_count, current_edges = parse_dimacs_edge(resume["edge_path"])
    if vertex_count != len(points):
        raise ValueError("resume vertex and certificate counts differ")
    old_coloring = tuple(int(color) for color in resume["whole_graph_coloring"])
    if not ColoringSAT(
        vertex_count,
        current_edges,
        5,
        break_color_symmetry=False,
    ).validate(old_coloring):
        raise ValueError("resume coloring is invalid")

    coordinates = np.asarray([point.approximate() for point in points])
    possible = cKDTree(coordinates).query_pairs(
        r=1.0 + tolerance,
        output_type="ndarray",
    )
    delta = coordinates[possible[:, 0]] - coordinates[possible[:, 1]]
    squared = np.einsum("ij,ij->i", delta, delta)
    annulus = possible[np.abs(squared - 1.0) <= 4.0 * tolerance]
    if progress is not None:
        progress(
            "numeric_annulus",
            {
                "pairs_within_radius": len(possible),
                "near_unit_pairs": len(annulus),
                "already_listed_edges": len(current_edges),
            },
        )
    added: set[tuple[int, int]] = set()
    for checked, (left_raw, right_raw) in enumerate(annulus, 1):
        left, right = int(left_raw), int(right_raw)
        edge = (left, right)
        if edge in current_edges:
            continue
        if exact.projective_unit(
            exact.coordinates[left],
            exact.coordinates[right],
        ):
            added.add(edge)
        if progress is not None and checked % 10_000 == 0:
            progress(
                "exact_interactions",
                {"checked": checked, "added": len(added)},
            )
    current_edges.update(added)
    checkpoint = _write_current_artifacts(
        output_directory=Path(output_directory),
        base_points=base_points,
        base_edges=base_edges,
        candidate_points=candidate_points,
        current_edges=current_edges,
        iterations=iterations,
        coloring=None,
        symmetry_clique=symmetry_clique,
        status="EXACT_INTERACTIONS_PENDING_SAT",
        complete_unit_distance_graph=False,
    )
    checkpoint["interaction_augmentation"] = {
        "source_report": str(resume_report_path),
        "numeric_pairs_within_one_plus_tolerance": len(possible),
        "numeric_near_unit_pairs": len(annulus),
        "new_exact_unit_edges": len(added),
        "numeric_filter_is_not_a_completeness_certificate": True,
    }
    checkpoint_report = (
        Path(str(checkpoint["graph_path"])).with_suffix("").with_suffix(".json")
    )
    _write_json(checkpoint_report, checkpoint)
    if progress is not None:
        progress(
            "checkpoint_written",
            {
                "vertices": vertex_count,
                "edges": len(current_edges),
                "new_exact_unit_edges": len(added),
                "cnf_path": checkpoint["cnf_path"],
            },
        )
    phases = [
        vertex * 5 + color + 1
        if color == old_coloring[vertex]
        else -(vertex * 5 + color + 1)
        for vertex in range(vertex_count)
        for color in range(5)
    ]
    solve_started = time.monotonic()
    coloring = ColoringSAT(
        vertex_count,
        current_edges,
        5,
        break_color_symmetry=False,
        symmetry_clique=symmetry_clique,
    ).solve(phases=phases)
    solve_seconds = time.monotonic() - solve_started
    status = (
        "EXACT_INTERACTION_SUBGRAPH_5_UNSAT_WITHOUT_PROOF"
        if coloring is None
        else "EXACT_INTERACTIONS_AUGMENTED_WITH_5_COLORING"
    )
    payload = _write_current_artifacts(
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
    for key in ("pair_cegis_rounds", "projective_closure_rounds", "duplicate_merge"):
        if key in resume:
            payload[key] = resume[key]
    payload["interaction_augmentation"] = {
        "source_report": str(resume_report_path),
        "numeric_pairs_within_one_plus_tolerance": len(possible),
        "numeric_near_unit_pairs": len(annulus),
        "new_exact_unit_edges": len(added),
        "solve_seconds": solve_seconds,
        "numeric_filter_is_not_a_completeness_certificate": True,
    }
    payload["elapsed_seconds"] = time.monotonic() - started
    report_path = Path(str(payload["graph_path"])).with_suffix("").with_suffix(".json")
    _write_json(report_path, payload)
    return payload
