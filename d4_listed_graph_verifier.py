"""Independent exact audit for listed-edge primary-D4 CEGIS artifacts.

The verifier imports the fixed algebra implementation but no blocker, pair, or
CEGIS search module.  It reconstructs every projective coordinate certificate,
checks every non-base edge as a 32-coefficient identity, validates the saved
coloring, and detects exact duplicate vertices after a numerical prefilter.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.spatial import cKDTree
import sympy as sp

from coloring_sat import ColoringSAT, parse_dimacs_edge
from d4_exact import D4ProjectivePoint, primary_d4_coordinate_backend
from exact_geometry import parse_mathematica_vertices
from graph_verifier import load_graph_json


def _coefficients(raw: object, key: str) -> tuple[sp.Rational, ...]:
    if not isinstance(raw, dict):
        raise ValueError("projective certificate must be an object")
    values = raw.get(key)
    if not isinstance(values, list) or len(values) != 32:
        raise ValueError(f"invalid coefficient vector {key}")
    return tuple(sp.Rational(str(value)) for value in values)


def _projective(raw: object) -> D4ProjectivePoint:
    return D4ProjectivePoint(
        _coefficients(raw, "numerator_x_coefficients"),
        _coefficients(raw, "numerator_y_coefficients"),
        _coefficients(raw, "denominator_coefficients"),
        -1,
    )


def verify_d4_listed_artifact(
    *,
    base_graph_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    report_path: str | Path,
) -> dict[str, object]:
    """Return exact geometry and coloring audit counts for one checkpoint."""

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
    with Path(report_path).open(encoding="utf-8") as handle:
        report = json.load(handle)
    candidates = [
        affine.algebra.normalize_projective(
            _projective(record["projective_certificate"])
        )
        for record in report["iterations"]
    ]
    projective = affine.as_projective_backend(candidates)
    vertex_count, edges = parse_dimacs_edge(report["edge_path"])
    if vertex_count != len(projective.coordinates):
        raise ValueError("edge file and coordinate certificate counts differ")
    if not base_edges.issubset(edges):
        raise ValueError("listed graph omits a base edge")

    nonbase_edges = sorted(edges - base_edges)
    invalid_edges = [
        (left, right)
        for left, right in nonbase_edges
        if not projective.projective_unit(
            projective.coordinates[left],
            projective.coordinates[right],
        )
    ]
    materialized = [
        *base_points,
        *(affine.materialize(point) for point in candidates),
    ]
    approximate = np.asarray(
        [point.approximate() for point in materialized],
        dtype=float,
    )
    if not np.isfinite(approximate).all():
        raise ValueError("non-finite coordinate approximation")
    possible_duplicates = cKDTree(approximate).query_pairs(
        r=1e-9,
        output_type="ndarray",
    )
    duplicate_pairs = [
        (int(left), int(right))
        for left, right in possible_duplicates
        if projective.projective_equal(
            projective.coordinates[int(left)],
            projective.coordinates[int(right)],
        )
    ]
    coloring_raw = report.get("whole_graph_coloring")
    coloring_valid = False
    class_sizes: list[int] | None = None
    if isinstance(coloring_raw, list):
        coloring = tuple(int(color) for color in coloring_raw)
        coloring_valid = ColoringSAT(
            vertex_count,
            edges,
            5,
            break_color_symmetry=False,
        ).validate(coloring)
        class_sizes = [coloring.count(color) for color in range(5)]
    return {
        "valid": not invalid_edges and not duplicate_pairs and coloring_valid,
        "vertices": vertex_count,
        "edges": len(edges),
        "base_edges": len(base_edges),
        "nonbase_edges_exactly_unit": len(nonbase_edges) - len(invalid_edges),
        "invalid_edges": invalid_edges,
        "numeric_near_duplicate_pairs": len(possible_duplicates),
        "exact_duplicate_pairs": duplicate_pairs,
        "coloring_valid": coloring_valid,
        "class_sizes": class_sizes,
    }
