"""Exact interactions among quadratic translations of the Heule graph."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Sequence

import sympy as sp

from coloring_sat import ColoringSAT, write_dimacs_edge
from exact_geometry import Point, parse_mathematica_vertices, serialize_graph
from graph_verifier import load_graph_json
from translated_copy_cegis import CircleTranslation, circle_translation


InteractionEdge = tuple[int, int, int, int]


def _canonical(value: sp.Expr | int) -> sp.Expr:
    return sp.radsimp(sp.cancel(sp.expand(value)))


@dataclass(frozen=True, slots=True)
class PairCertificate:
    left_copy: int
    right_copy: int
    method: str
    edge_count: int
    complete: bool
    same_extension: bool
    parallel_candidate_pairs: int
    separating_coefficient: str


@dataclass(frozen=True, slots=True)
class InteractionResult:
    edges: frozenset[InteractionEdge]
    pair_certificates: tuple[PairCertificate, ...]


@dataclass(frozen=True, slots=True)
class StageArtifact:
    stage: int
    edge_count: int
    added_interaction_edges: int
    graph_path: Path
    edge_path: Path
    cnf_path: Path
    report_path: Path


def _translated(point: Point, translation: CircleTranslation) -> Point:
    return Point(
        point.x + translation.point.x,
        point.y + translation.point.y,
    )


def _same_extension_edges(
    source_points: Sequence[Point],
    left_copy: int,
    right_copy: int,
    left_r: Point,
    right_r: Point,
    left_translation: CircleTranslation,
    right_translation: CircleTranslation,
) -> tuple[set[InteractionEdge], int]:
    """Enumerate one pair by separating the coefficients of ``1,h``."""

    s = Point(right_r.x - left_r.x, right_r.y - left_r.y)
    s_squared = s.squared_distance(Point(0, 0))
    if s_squared == 0:
        raise ValueError("distinct copies use the same translation center")

    parallel_lines: dict[sp.Expr, list[int]] = defaultdict(list)
    for index, point in enumerate(source_points):
        cross_coordinate = _canonical(point.x * s.y - point.y * s.x)
        parallel_lines[cross_coordinate].append(index)

    edges: set[InteractionEdge] = set()
    parallel_candidate_pairs = 0
    h_squared = left_translation.h_squared
    for indices in parallel_lines.values():
        parallel_candidate_pairs += len(indices) * len(indices)
        for left_index in indices:
            left_point = source_points[left_index]
            for right_index in indices:
                right_point = source_points[right_index]
                wx = right_point.x - left_point.x
                wy = right_point.y - left_point.y
                constant_coefficient = _canonical(
                    wx * wx
                    + wy * wy
                    + wx * s.x
                    + wy * s.y
                    + s_squared * (sp.Rational(1, 4) + h_squared)
                    - 1
                )
                if constant_coefficient != 0:
                    continue
                translated_left = _translated(left_point, left_translation)
                translated_right = _translated(right_point, right_translation)
                if translated_left.squared_distance(translated_right) != 1:
                    raise AssertionError(
                        "basis-separated same-extension edge is not unit"
                    )
                edges.add(
                    (left_copy, left_index, right_copy, right_index)
                )
    return edges, parallel_candidate_pairs


def complete_copy_interactions(
    source_points: Sequence[Point],
    selected_r_indices: Sequence[int],
) -> InteractionResult:
    """Return every exact unit edge among all selected translated copies.

    Same-extension pairs are enumerated by exact parallel-line grouping.
    Mixed-extension pairs are excluded by their nonzero ``h_+ h_-``
    coefficient.
    """

    if len(set(selected_r_indices)) != len(selected_r_indices):
        raise ValueError("selected translation indices must be distinct")
    translations = [
        circle_translation(source_points[index], branch=1)
        for index in selected_r_indices
    ]
    all_edges: set[InteractionEdge] = set()
    certificates: list[PairCertificate] = []

    for left_copy in range(len(selected_r_indices)):
        left_r = source_points[selected_r_indices[left_copy]]
        left_translation = translations[left_copy]
        for right_copy in range(left_copy + 1, len(selected_r_indices)):
            right_r = source_points[selected_r_indices[right_copy]]
            right_translation = translations[right_copy]
            same_extension = (
                _canonical(
                    left_translation.h_squared
                    - right_translation.h_squared
                )
                == 0
            )
            if same_extension:
                pair_edges, candidate_count = _same_extension_edges(
                    source_points,
                    left_copy,
                    right_copy,
                    left_r,
                    right_r,
                    left_translation,
                    right_translation,
                )
                all_edges.update(pair_edges)
                certificates.append(
                    PairCertificate(
                        left_copy=left_copy,
                        right_copy=right_copy,
                        method="same_extension_basis_separation",
                        edge_count=len(pair_edges),
                        complete=True,
                        same_extension=True,
                        parallel_candidate_pairs=candidate_count,
                        separating_coefficient="0 (parallel-line grouping)",
                    )
                )
                continue

            product = _canonical(
                left_translation.h_squared
                * right_translation.h_squared
            )
            if product != sp.Rational(41, 8):
                raise ValueError(
                    "mixed extensions do not have the certified product 41/8"
                )
            dot_product = _canonical(
                left_r.x * right_r.x + left_r.y * right_r.y
            )
            if dot_product == 0:
                determinant = _canonical(
                    left_r.x * right_r.y - left_r.y * right_r.x
                )
                if determinant == 0:
                    raise ValueError(
                        "mixed-extension vectors are parallel with zero dot "
                        "product, contradicting their nonzero lengths"
                    )
                left_squared = left_r.squared_distance(Point(0, 0))
                right_squared = right_r.squared_distance(Point(0, 0))
                residual_constant = _canonical(
                    left_squared * left_translation.h_squared
                    + right_squared * right_translation.h_squared
                    - 1
                )
                if residual_constant == 0:
                    raise ValueError(
                        "orthogonal mixed-extension pair has a vanishing "
                        "residual constant and needs source-point enumeration"
                    )
                separating_coefficient = (
                    "h_product=0; orthogonal residual="
                    + sp.sstr(residual_constant)
                )
            else:
                separating_coefficient = sp.sstr(-2 * dot_product)
            certificates.append(
                PairCertificate(
                    left_copy=left_copy,
                    right_copy=right_copy,
                    method="mixed_extension_product_coefficient",
                    edge_count=0,
                    complete=True,
                    same_extension=False,
                    parallel_candidate_pairs=0,
                    separating_coefficient=separating_coefficient,
                )
            )

    return InteractionResult(
        edges=frozenset(all_edges),
        pair_certificates=tuple(certificates),
    )


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def export_interaction_stages(
    *,
    graph_path: str | Path,
    report_path: str | Path,
    source_vertex_path: str | Path,
    output_directory: str | Path,
    color_count: int = 5,
    symmetry_clique: Sequence[int] = (548, 1149, 668),
) -> list[StageArtifact]:
    """Export cumulative nonempty interaction stages as graph/EDGE/CNF."""

    graph_path = Path(graph_path)
    report_path = Path(report_path)
    output_directory = Path(output_directory)
    points, base_edges = load_graph_json(graph_path)
    with report_path.open(encoding="utf-8") as handle:
        source_report = json.load(handle)
    selected = [
        int(index)
        for index in source_report["selected_source_r_indices_zero_based"]
    ]
    source_points = parse_mathematica_vertices(source_vertex_path)
    source_count = len(source_points)
    base_vertex_count = len(points) - len(selected) * source_count
    if base_vertex_count <= 0:
        raise ValueError("graph size is inconsistent with translated copies")
    if len(points) != 8329 or len(base_edges) != 51378:
        raise ValueError("unexpected primary-D4 source artifact")

    interactions = complete_copy_interactions(source_points, selected)
    edges_by_pair: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
    for left_copy, left_vertex, right_copy, right_vertex in interactions.edges:
        left = base_vertex_count + left_copy * source_count + left_vertex
        right = base_vertex_count + right_copy * source_count + right_vertex
        edges_by_pair[(left_copy, right_copy)].add((left, right))
    nonempty_pairs = sorted(edges_by_pair)
    if [len(edges_by_pair[pair]) for pair in nonempty_pairs] != [648, 655, 651]:
        raise AssertionError("unexpected nonempty interaction pair counts")

    output_directory.mkdir(parents=True, exist_ok=True)
    cumulative_edges = set(base_edges)
    artifacts: list[StageArtifact] = []
    cumulative_interactions = 0
    for stage, pair in enumerate(nonempty_pairs, 1):
        pair_edges = edges_by_pair[pair]
        cumulative_edges.update(pair_edges)
        cumulative_interactions += len(pair_edges)
        stem = f"heule8329-interactions-stage{stage}"
        stage_graph = output_directory / f"{stem}.graph.json"
        stage_edge = output_directory / f"{stem}.edge"
        stage_cnf = output_directory / f"{stem}-5.cnf"
        stage_report = output_directory / f"{stem}.json"

        _write_json(
            stage_graph,
            serialize_graph(
                f"8329-point translated-copy interaction stage {stage}",
                points,
                cumulative_edges,
            ),
        )
        write_dimacs_edge(stage_edge, len(points), cumulative_edges)
        encoding = ColoringSAT(
            len(points),
            cumulative_edges,
            color_count,
            break_color_symmetry=False,
            symmetry_clique=symmetry_clique,
        )
        encoding.write_dimacs_cnf(stage_cnf)
        payload: dict[str, object] = {
            "status": "UNSOLVED",
            "lower_bound_6_proved": False,
            "stage": stage,
            "vertices": len(points),
            "edges": len(cumulative_edges),
            "added_interaction_edges": cumulative_interactions,
            "new_pair": list(pair),
            "new_pair_edges": len(pair_edges),
            "complete_unit_distance_graph": stage == len(nonempty_pairs),
            "symmetry_clique_zero_based": list(symmetry_clique),
            "cnf_variables": encoding.variable_count,
            "cnf_clauses": len(encoding.clauses),
            "graph_path": str(stage_graph),
            "edge_path": str(stage_edge),
            "cnf_path": str(stage_cnf),
            "source_graph_path": str(graph_path),
            "source_report_path": str(report_path),
            "pair_certificates": [
                {
                    "left_copy": certificate.left_copy,
                    "right_copy": certificate.right_copy,
                    "method": certificate.method,
                    "edge_count": certificate.edge_count,
                    "complete": certificate.complete,
                    "same_extension": certificate.same_extension,
                    "parallel_candidate_pairs": (
                        certificate.parallel_candidate_pairs
                    ),
                    "separating_coefficient": (
                        certificate.separating_coefficient
                    ),
                }
                for certificate in interactions.pair_certificates
            ],
        }
        _write_json(stage_report, payload)
        payload["sha256"] = {
            "graph": _sha256(stage_graph),
            "edge": _sha256(stage_edge),
            "cnf": _sha256(stage_cnf),
        }
        _write_json(stage_report, payload)
        artifacts.append(
            StageArtifact(
                stage=stage,
                edge_count=len(cumulative_edges),
                added_interaction_edges=cumulative_interactions,
                graph_path=stage_graph,
                edge_path=stage_edge,
                cnf_path=stage_cnf,
                report_path=stage_report,
            )
        )
    return artifacts
