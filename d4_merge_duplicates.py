"""Merge exactly coincident projective vertices in a D4 CEGIS checkpoint."""

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


_SCALAR_INDEX_KEYS = {
    "candidate_vertex_zero_based",
    "generator_4_zero_based",
    "fifth_neighbor_zero_based",
}
_LIST_INDEX_KEYS = {
    "generator_01_zero_based",
    "generator_23_zero_based",
    "listed_base_neighbors_zero_based",
    "listed_neighbors_zero_based",
    "complete_base_neighbors_zero_based",
    "candidate_interaction_neighbors_zero_based",
}


def merge_exact_duplicates(
    *,
    base_graph_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    resume_report_path: str | Path,
    output_directory: str | Path,
    symmetry_clique: Sequence[int] = (548, 1149, 668),
    progress: Callable[[str, dict[str, object]], None] | None = None,
) -> dict[str, object]:
    """Identify projective equalities, quotient the graph, and re-solve it."""

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
    projective = affine.as_projective_backend(candidates)
    old_vertex_count, old_edges = parse_dimacs_edge(resume["edge_path"])
    if old_vertex_count != len(points):
        raise ValueError("resume vertex and coordinate counts differ")

    approximate = np.asarray([point.approximate() for point in points])
    possible = cKDTree(approximate).query_pairs(r=1e-9, output_type="ndarray")
    parent = list(range(old_vertex_count))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    duplicate_pairs: list[tuple[int, int]] = []
    for left_raw, right_raw in possible:
        left, right = int(left_raw), int(right_raw)
        if not projective.projective_equal(
            projective.coordinates[left],
            projective.coordinates[right],
        ):
            continue
        if left < len(base_points) or right < len(base_points):
            raise ValueError("candidate duplicates a base vertex")
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root
            duplicate_pairs.append((left, right))
    if not duplicate_pairs:
        raise ValueError("checkpoint has no exact duplicates")

    representative_to_new: dict[int, int] = {}
    old_to_new: list[int] = []
    for old_vertex in range(old_vertex_count):
        representative = find(old_vertex)
        if representative not in representative_to_new:
            representative_to_new[representative] = len(representative_to_new)
        old_to_new.append(representative_to_new[representative])
    if old_to_new[: len(base_points)] != list(range(len(base_points))):
        raise AssertionError("base vertex numbering changed")
    new_edges = {
        (min(old_to_new[left], old_to_new[right]),
         max(old_to_new[left], old_to_new[right]))
        for left, right in old_edges
        if old_to_new[left] != old_to_new[right]
    }

    kept_candidate_indices = [
        index
        for index in range(len(candidates))
        if find(len(base_points) + index) == len(base_points) + index
    ]
    new_candidates = [candidates[index] for index in kept_candidate_indices]
    new_candidate_points = [candidate_points[index] for index in kept_candidate_indices]

    def remap_record(record: dict[str, object]) -> dict[str, object]:
        remapped = json.loads(json.dumps(record))
        for key in _SCALAR_INDEX_KEYS:
            if key in remapped and isinstance(remapped[key], int) and remapped[key] >= 0:
                remapped[key] = old_to_new[remapped[key]]
        for key in _LIST_INDEX_KEYS:
            if key in remapped and isinstance(remapped[key], list):
                remapped[key] = sorted(
                    {old_to_new[int(vertex)] for vertex in remapped[key]}
                )
        certificate = remapped.get("projective_certificate")
        if isinstance(certificate, dict):
            fifth = certificate.get("fifth_neighbor_zero_based")
            if isinstance(fifth, int) and fifth >= 0:
                certificate["fifth_neighbor_zero_based"] = old_to_new[fifth]
        return remapped

    new_iterations = [
        remap_record(iterations[index]) for index in kept_candidate_indices
    ]
    for iteration_number, record in enumerate(new_iterations, 1):
        record["iteration"] = iteration_number
    old_coloring = tuple(int(color) for color in resume["whole_graph_coloring"])
    phase_coloring = [0] * len(representative_to_new)
    for old_vertex, new_vertex in enumerate(old_to_new):
        if old_vertex == find(old_vertex):
            phase_coloring[new_vertex] = old_coloring[old_vertex]
    phases = [
        vertex * 5 + color + 1
        if color == phase_coloring[vertex]
        else -(vertex * 5 + color + 1)
        for vertex in range(len(phase_coloring))
        for color in range(5)
    ]
    if progress is not None:
        progress(
            "duplicates_merged",
            {
                "old_vertices": old_vertex_count,
                "new_vertices": len(representative_to_new),
                "duplicate_pairs": duplicate_pairs,
                "old_edges": len(old_edges),
                "new_edges": len(new_edges),
            },
        )
    coloring = ColoringSAT(
        len(representative_to_new),
        new_edges,
        5,
        break_color_symmetry=False,
        symmetry_clique=symmetry_clique,
    ).solve(phases=phases)
    status = (
        "EXACT_DUPLICATE_QUOTIENT_5_UNSAT_WITHOUT_PROOF"
        if coloring is None
        else "EXACT_DUPLICATES_MERGED_WITH_5_COLORING"
    )
    payload = _write_current_artifacts(
        output_directory=Path(output_directory),
        base_points=base_points,
        base_edges=base_edges,
        candidate_points=new_candidate_points,
        current_edges=new_edges,
        iterations=new_iterations,
        coloring=coloring,
        symmetry_clique=symmetry_clique,
        status=status,
        complete_unit_distance_graph=False,
    )
    for key in ("pair_cegis_rounds", "projective_closure_rounds"):
        if key in resume:
            payload[key] = resume[key]
    payload["duplicate_merge"] = {
        "source_report": str(resume_report_path),
        "duplicate_pairs_zero_based": [list(pair) for pair in duplicate_pairs],
        "removed_vertices": old_vertex_count - len(representative_to_new),
        "removed_or_coalesced_edges": len(old_edges) - len(new_edges),
    }
    payload["elapsed_seconds"] = time.monotonic() - started
    report_path = Path(str(payload["graph_path"])).with_suffix("").with_suffix(".json")
    _write_json(report_path, payload)
    return payload
