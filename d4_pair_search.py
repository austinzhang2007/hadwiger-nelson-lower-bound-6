"""Exact search for adjacent four-color-forced candidate pairs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.spatial import cKDTree

from blocking_point_search import (
    _cross_color_center_pairs,
    _partition_intersections,
    _structured_rows,
)
from d4_exact import D4CoordinateBackend, D4ProjectivePoint
from exact_geometry import Point


@dataclass(frozen=True, slots=True)
class ForcedPoint:
    point: Point
    projective: D4ProjectivePoint
    neighbors: tuple[int, ...]
    forced_color: int
    generator_pairs: tuple[tuple[int, int], tuple[int, int]]


@dataclass(frozen=True, slots=True)
class ForcedPointStats:
    missing_color: int
    center_pairs: int
    first_partition_pairs: int
    second_partition_pairs: int
    common_intersections: int
    nonexisting_intersections: int
    exact_candidates: int


def exact_four_color_forced_points(
    points: Sequence[Point],
    coloring: Sequence[int],
    backend: D4CoordinateBackend,
    *,
    missing_color: int,
    quantization_digits: int = 10,
    tolerance: float = 1e-8,
) -> tuple[list[ForcedPoint], ForcedPointStats]:
    """Return exact points with listed neighbors in all colors but one."""

    if len(points) != len(coloring) or len(backend.coordinates) != len(points):
        raise ValueError("point, coloring, and exact-backend lengths differ")
    if missing_color not in range(5):
        raise ValueError("missing_color must be in range(5)")
    present = [color for color in range(5) if color != missing_color]
    coordinates = np.asarray(
        [point.approximate() for point in points],
        dtype=float,
    )
    color_array = np.asarray(coloring, dtype=np.int8)
    first_center_pairs = _cross_color_center_pairs(
        coordinates,
        color_array,
        present[0],
        present[1],
        2.0 + tolerance,
    )
    second_center_pairs = _cross_color_center_pairs(
        coordinates,
        color_array,
        present[2],
        present[3],
        2.0 + tolerance,
    )
    scale = 10**quantization_digits
    first_keys, first_points, first_generators, first_pair_count = (
        _partition_intersections(
            coordinates,
            color_array,
            first_center_pairs,
            present[0],
            present[1],
            scale,
            tolerance,
        )
    )
    second_keys, _, second_generators, second_pair_count = (
        _partition_intersections(
            coordinates,
            color_array,
            second_center_pairs,
            present[2],
            present[3],
            scale,
            tolerance,
        )
    )
    _, first_indices, second_indices = np.intersect1d(
        _structured_rows(first_keys),
        _structured_rows(second_keys),
        assume_unique=True,
        return_indices=True,
    )
    common_count = len(first_indices)
    existing_keys = np.rint(coordinates * scale).astype(np.int64)
    existing = _structured_rows(np.unique(existing_keys, axis=0))
    is_existing = np.isin(
        _structured_rows(first_keys[first_indices]),
        existing,
        assume_unique=False,
    )
    first_indices = first_indices[~is_existing]
    second_indices = second_indices[~is_existing]

    recovered: dict[tuple[int, int], ForcedPoint] = {}
    for first_index, second_index in zip(first_indices, second_indices):
        first_generator = tuple(
            int(index) for index in first_generators[first_index]
        )
        second_generator = tuple(
            int(index) for index in second_generators[second_index]
        )
        projective = backend.recover_four_projective(
            first_generator,
            second_generator,
        )
        if projective is None:
            continue
        point = backend.materialize(projective)
        approximate = point.approximate()
        numeric_key = (
            round(approximate[0], quantization_digits),
            round(approximate[1], quantization_digits),
        )
        neighbors = tuple(sorted({*first_generator, *second_generator}))
        if {
            int(coloring[index]) for index in neighbors
        } != set(present):
            continue
        recovered.setdefault(
            numeric_key,
            ForcedPoint(
                point=point,
                projective=projective,
                neighbors=neighbors,
                forced_color=missing_color,
                generator_pairs=(first_generator, second_generator),
            ),
        )
    candidates = sorted(
        recovered.values(),
        key=lambda item: item.point.approximate(),
    )
    return candidates, ForcedPointStats(
        missing_color=missing_color,
        center_pairs=len(first_center_pairs) + len(second_center_pairs),
        first_partition_pairs=first_pair_count,
        second_partition_pairs=second_pair_count,
        common_intersections=common_count,
        nonexisting_intersections=len(first_indices),
        exact_candidates=len(candidates),
    )


def exact_forced_pair_edges(
    candidates: Sequence[ForcedPoint],
    backend: D4CoordinateBackend,
    *,
    tolerance: float = 1e-8,
) -> set[tuple[int, int]]:
    """Return every exact same-forced-color unit pair after numeric filtering."""

    if len(candidates) < 2:
        return set()
    coordinates = np.asarray(
        [candidate.point.approximate() for candidate in candidates],
        dtype=float,
    )
    possible = cKDTree(coordinates).query_pairs(
        r=1.0 + tolerance,
        output_type="ndarray",
    )
    edges: set[tuple[int, int]] = set()
    for left, right in possible:
        left_index, right_index = int(left), int(right)
        if candidates[left_index].forced_color != candidates[right_index].forced_color:
            continue
        squared = float(np.sum((coordinates[left_index] - coordinates[right_index]) ** 2))
        if abs(squared - 1.0) > 4.0 * tolerance:
            continue
        if backend.projective_unit(
            candidates[left_index].projective,
            candidates[right_index].projective,
        ):
            edges.add((left_index, right_index))
    return edges
