"""Independent verifier for the complete translated-copy interaction graph.

This module deliberately does not import either search implementation.
"""

from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
from typing import Sequence

import sympy as sp

from exact_geometry import Point, parse_mathematica_vertices
from graph_verifier import load_graph_json
from translated_copy_verifier import _circle_translation, _parse_edge_file


def _canonical(value: sp.Expr | int) -> sp.Expr:
    return sp.radsimp(sp.cancel(sp.expand(value)))


def _h_squared(r: Point) -> sp.Expr:
    squared = r.squared_distance(Point(0, 0))
    return _canonical((4 - squared) / (4 * squared))


def _same_extension_interactions(
    source: Sequence[Point],
    left_copy: int,
    right_copy: int,
    left_r: Point,
    right_r: Point,
    left_translation: Point,
    right_translation: Point,
    h_squared: sp.Expr,
) -> set[tuple[int, int]]:
    s = Point(right_r.x - left_r.x, right_r.y - left_r.y)
    s_squared = s.squared_distance(Point(0, 0))
    lines: dict[sp.Expr, list[int]] = defaultdict(list)
    for index, point in enumerate(source):
        lines[_canonical(point.x * s.y - point.y * s.x)].append(index)

    result: set[tuple[int, int]] = set()
    for indices in lines.values():
        for left_index in indices:
            left = source[left_index]
            for right_index in indices:
                right = source[right_index]
                wx, wy = right.x - left.x, right.y - left.y
                constant = _canonical(
                    wx * wx
                    + wy * wy
                    + wx * s.x
                    + wy * s.y
                    + s_squared * (sp.Rational(1, 4) + h_squared)
                    - 1
                )
                if constant != 0:
                    continue
                shifted_left = Point(
                    left.x + left_translation.x,
                    left.y + left_translation.y,
                )
                shifted_right = Point(
                    right.x + right_translation.x,
                    right.y + right_translation.y,
                )
                if shifted_left.squared_distance(shifted_right) != 1:
                    raise AssertionError(
                        f"copies {left_copy},{right_copy}: non-unit edge"
                    )
                result.add((left_index, right_index))
    return result


def verify_complete_interaction_artifact(
    *,
    graph_path: str | Path,
    stage_report_path: str | Path,
    primary_report_path: str | Path,
    base_graph_path: str | Path,
    source_vertex_path: str | Path,
    source_edge_path: str | Path,
) -> dict[str, object]:
    """Reconstruct all 8329 coordinates and every unit edge independently."""

    points, listed_edges = load_graph_json(graph_path)
    base_points, base_edges = load_graph_json(base_graph_path)
    source = parse_mathematica_vertices(source_vertex_path)
    source_count, source_edges = _parse_edge_file(source_edge_path)
    with Path(primary_report_path).open(encoding="utf-8") as handle:
        primary_report = json.load(handle)
    with Path(stage_report_path).open(encoding="utf-8") as handle:
        stage_report = json.load(handle)
    selected = [
        int(index)
        for index in primary_report["selected_source_r_indices_zero_based"]
    ]
    if source_count != len(source):
        raise ValueError("source coordinate and edge counts disagree")

    translations: list[Point] = []
    h_squares: list[sp.Expr] = []
    expected_points = list(base_points)
    translation_certificates_valid = True
    for source_r_index in selected:
        r = source[source_r_index]
        translation, _, nonabelian = _circle_translation(r)
        translations.append(translation)
        h_squares.append(_h_squared(r))
        translation_certificates_valid &= bool(
            nonabelian
            and translation.squared_distance(Point(0, 0)) == 1
            and translation.squared_distance(r) == 1
        )
        expected_points.extend(
            Point(point.x + translation.x, point.y + translation.y)
            for point in source
        )

    coordinates_match = points == expected_points
    duplicate_vertices = len(points) - len(set(points))
    expected_edges = set(base_edges)
    base_lookup = {point: index for index, point in enumerate(base_points)}
    for copy, source_r_index in enumerate(selected):
        offset = len(base_points) + copy * source_count
        r = source[source_r_index]
        expected_edges.update(
            (offset + left, offset + right) for left, right in source_edges
        )
        for source_index, point in enumerate(source):
            for anchor in (
                point,
                Point(point.x + r.x, point.y + r.y),
            ):
                base_index = base_lookup.get(anchor)
                if base_index is not None:
                    expected_edges.add((base_index, offset + source_index))

    interaction_edges = 0
    pair_certificates = 0
    mixed_orthogonal_pairs = 0
    for left_copy in range(len(selected)):
        left_r = source[selected[left_copy]]
        for right_copy in range(left_copy + 1, len(selected)):
            pair_certificates += 1
            right_r = source[selected[right_copy]]
            if _canonical(h_squares[left_copy] - h_squares[right_copy]) == 0:
                pair_edges = _same_extension_interactions(
                    source,
                    left_copy,
                    right_copy,
                    left_r,
                    right_r,
                    translations[left_copy],
                    translations[right_copy],
                    h_squares[left_copy],
                )
                left_offset = len(base_points) + left_copy * source_count
                right_offset = len(base_points) + right_copy * source_count
                expected_edges.update(
                    (left_offset + left, right_offset + right)
                    for left, right in pair_edges
                )
                interaction_edges += len(pair_edges)
                continue

            if _canonical(
                h_squares[left_copy] * h_squares[right_copy]
            ) != sp.Rational(41, 8):
                raise AssertionError("unexpected mixed-extension product")
            dot = _canonical(
                left_r.x * right_r.x + left_r.y * right_r.y
            )
            if dot != 0:
                continue
            determinant = _canonical(
                left_r.x * right_r.y - left_r.y * right_r.x
            )
            if determinant == 0:
                raise AssertionError("zero mixed coefficients are dependent")
            residual = _canonical(
                left_r.squared_distance(Point(0, 0))
                * h_squares[left_copy]
                + right_r.squared_distance(Point(0, 0))
                * h_squares[right_copy]
                - 1
            )
            if residual == 0:
                raise AssertionError("orthogonal mixed residual vanished")
            mixed_orthogonal_pairs += 1

    base_edges_are_unit = all(
        base_points[left].squared_distance(base_points[right]) == 1
        for left, right in base_edges
    )
    source_edges_are_unit = all(
        source[left].squared_distance(source[right]) == 1
        for left, right in source_edges
    )
    coloring_raw = stage_report.get("coloring")
    coloring_conflicts = None
    if coloring_raw is not None:
        coloring = tuple(int(color) for color in coloring_raw)
        coloring_conflicts = (
            None
            if len(coloring) != len(points)
            else sum(
                coloring[left] == coloring[right]
                for left, right in listed_edges
            )
        )
    missing = expected_edges - listed_edges
    unexpected = listed_edges - expected_edges
    valid = bool(
        coordinates_match
        and duplicate_vertices == 0
        and translation_certificates_valid
        and base_edges_are_unit
        and source_edges_are_unit
        and not missing
        and not unexpected
        and (coloring_conflicts is None or coloring_conflicts == 0)
        and stage_report.get("complete_unit_distance_graph") is True
    )
    return {
        "valid": valid,
        "vertices": len(points),
        "edges": len(listed_edges),
        "interaction_edges": interaction_edges,
        "pair_certificates": pair_certificates,
        "mixed_orthogonal_pairs": mixed_orthogonal_pairs,
        "coordinates_match": coordinates_match,
        "duplicate_vertices": duplicate_vertices,
        "translation_certificates_valid": translation_certificates_valid,
        "base_edges_symbolically_unit": base_edges_are_unit,
        "source_edges_symbolically_unit": source_edges_are_unit,
        "missing_edges": len(missing),
        "unexpected_edges": len(unexpected),
        "coloring_conflicts": coloring_conflicts,
    }
