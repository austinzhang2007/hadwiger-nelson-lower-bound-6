"""Color-partitioned numerical search with exact blocker recovery.

Numerical coordinates are used only to find seeds.  Every returned candidate
and every recorded unit edge is reconstructed and checked symbolically.
"""

from __future__ import annotations

from dataclasses import dataclass
import time
from typing import TYPE_CHECKING, Callable, Sequence

import numpy as np
from scipy.spatial import cKDTree

from exact_geometry import Point, _field_coordinates

if TYPE_CHECKING:
    from d4_exact import D4CoordinateBackend


@dataclass(frozen=True, slots=True)
class ExactBlockingPoint:
    point: Point
    neighbors: tuple[int, ...]
    neighbor_colors: tuple[int, ...]
    generator_pair: tuple[int, int]


@dataclass(frozen=True, slots=True)
class BlockingSearchStats:
    center_pairs: int
    first_partition_pairs: int
    second_partition_pairs: int
    first_distinct_intersections: int
    second_distinct_intersections: int
    common_partition_intersections: int
    nonexisting_common_intersections: int
    fifth_color_numeric_hits: int
    exact_blockers: int
    elapsed_seconds: float


def _structured_rows(values: np.ndarray) -> np.ndarray:
    contiguous = np.ascontiguousarray(values, dtype=np.int64)
    return contiguous.view([("x", np.int64), ("y", np.int64)]).reshape(-1)


def _cross_color_center_pairs(
    coordinates: np.ndarray,
    color_array: np.ndarray,
    left_color: int,
    right_color: int,
    radius: float,
) -> np.ndarray:
    """Return only cross-color center pairs within ``radius``.

    Constructing every nearby pair and discarding roughly 90% by color became
    the dominant memory cost after projective closure exceeded 20k vertices.
    The sparse bipartite query is equivalent for a fixed color partition.
    """

    left_indices = np.flatnonzero(color_array == left_color)
    right_indices = np.flatnonzero(color_array == right_color)
    if not len(left_indices) or not len(right_indices):
        return np.empty((0, 2), dtype=np.int64)
    distances = cKDTree(coordinates[left_indices]).sparse_distance_matrix(
        cKDTree(coordinates[right_indices]),
        max_distance=radius,
        output_type="coo_matrix",
    )
    return np.column_stack(
        (left_indices[distances.row], right_indices[distances.col])
    ).astype(np.int64, copy=False)


def _partition_intersections(
    coordinates: np.ndarray,
    color_array: np.ndarray,
    center_pairs: np.ndarray,
    left_color: int,
    right_color: int,
    scale: int,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    left = center_pairs[:, 0]
    right = center_pairs[:, 1]
    mask = (
        ((color_array[left] == left_color) & (color_array[right] == right_color))
        | ((color_array[left] == right_color) & (color_array[right] == left_color))
    )
    selected = center_pairs[mask]
    if not len(selected):
        empty_keys = np.empty((0, 2), dtype=np.int64)
        empty_points = np.empty((0, 2), dtype=float)
        return empty_keys, empty_points, selected, 0

    delta = coordinates[selected[:, 1]] - coordinates[selected[:, 0]]
    distance_squared = np.einsum("ij,ij->i", delta, delta)
    valid = (
        (distance_squared > tolerance * tolerance)
        & (distance_squared <= 4.0 + tolerance)
    )
    selected = selected[valid]
    delta = delta[valid]
    distance_squared = distance_squared[valid]
    midpoint = (
        coordinates[selected[:, 0]] + coordinates[selected[:, 1]]
    ) / 2.0
    height_over_distance = np.sqrt(
        np.maximum(0.0, 1.0 - distance_squared / 4.0) / distance_squared
    )
    perpendicular = np.column_stack((-delta[:, 1], delta[:, 0]))
    displacement = perpendicular * height_over_distance[:, None]
    raw_points = np.concatenate(
        (midpoint - displacement, midpoint + displacement),
        axis=0,
    )
    raw_generators = np.concatenate((selected, selected), axis=0)
    raw_keys = np.rint(raw_points * scale).astype(np.int64)
    keys, first_indices = np.unique(raw_keys, axis=0, return_index=True)
    return (
        keys,
        raw_points[first_indices],
        raw_generators[first_indices],
        len(selected),
    )


def exact_five_color_blockers(
    points: Sequence[Point],
    coloring: Sequence[int],
    *,
    quantization_digits: int = 10,
    tolerance: float = 1e-8,
    progress: Callable[[str, dict[str, int]], None] | None = None,
    exact_backend: "D4CoordinateBackend | None" = None,
    max_blockers: int | None = None,
) -> tuple[list[ExactBlockingPoint], BlockingSearchStats]:
    """Find candidates adjacent to all five colors and certify them exactly.

    If a point sees all five colors, it is simultaneously an intersection of a
    color-0/color-1 circle pair and of a color-2/color-3 circle pair.  This
    necessary condition makes the numerical search much smaller while retaining
    every single-point blocker for the supplied coloring.
    """

    started = time.monotonic()
    if len(points) != len(coloring):
        raise ValueError("point and coloring lengths differ")
    if max_blockers is not None and max_blockers < 1:
        raise ValueError("max_blockers must be positive")
    if set(int(color) for color in coloring) != set(range(5)):
        raise ValueError("the search requires a surjective five-coloring")
    coordinates = np.asarray(
        [point.approximate() for point in points],
        dtype=float,
    )
    color_array = np.asarray(coloring, dtype=np.int8)
    first_center_pairs = _cross_color_center_pairs(
        coordinates,
        color_array,
        0,
        1,
        2.0 + tolerance,
    )
    second_center_pairs = _cross_color_center_pairs(
        coordinates,
        color_array,
        2,
        3,
        2.0 + tolerance,
    )
    scale = 10**quantization_digits
    first_keys, first_points, first_generators, first_pair_count = (
        _partition_intersections(
            coordinates,
            color_array,
            first_center_pairs,
            0,
            1,
            scale,
            tolerance,
        )
    )
    second_keys, second_points, second_generators, second_pair_count = (
        _partition_intersections(
            coordinates,
            color_array,
            second_center_pairs,
        2,
        3,
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
    is_existing = np.isin(
        _structured_rows(first_keys[first_indices]),
        _structured_rows(np.unique(existing_keys, axis=0)),
        assume_unique=False,
    )
    first_indices = first_indices[~is_existing]
    second_indices = second_indices[~is_existing]
    if progress is not None:
        progress(
            "numeric_intersection",
            {
                "center_pairs": len(first_center_pairs) + len(second_center_pairs),
                "first_partition_pairs": first_pair_count,
                "second_partition_pairs": second_pair_count,
                "first_distinct_intersections": len(first_keys),
                "second_distinct_intersections": len(second_keys),
                "common_partition_intersections": common_count,
                "nonexisting_common_intersections": len(first_indices),
            },
        )

    fifth_indices = np.flatnonzero(color_array == 4)
    fifth_coordinates = coordinates[fifth_indices]
    fifth_tree = cKDTree(fifth_coordinates)
    numeric_hits: list[tuple[int, int, tuple[int, ...]]] = []
    for first_index, second_index in zip(first_indices, second_indices):
        approximate = first_points[first_index]
        possible = fifth_tree.query_ball_point(
            approximate,
            r=1.0 + tolerance,
        )
        near_unit = tuple(
            int(fifth_indices[local_index])
            for local_index in possible
            if abs(
                float(
                    np.sum(
                        (
                            fifth_coordinates[local_index]
                            - approximate
                        )
                        ** 2
                    )
                )
                - 1.0
            )
            <= 4.0 * tolerance
        )
        if near_unit:
            numeric_hits.append(
                (int(first_index), int(second_index), near_unit)
            )
    if progress is not None:
        progress(
            "fifth_color_filter",
            {"fifth_color_numeric_hits": len(numeric_hits)},
        )

    if exact_backend is not None and len(exact_backend.coordinates) != len(points):
        raise ValueError("exact backend and point counts differ")
    existing = set(points)
    recovered: dict[Point, ExactBlockingPoint] = {}
    if exact_backend is not None:
        certified: list[
            tuple[
                tuple[int, int, int],
                tuple[int, int],
                tuple[int, int],
                object,
            ]
        ] = []
        for first_index, second_index, fifth_neighbors in numeric_hits:
            first_generator = tuple(
                int(index) for index in first_generators[first_index]
            )
            second_generator = tuple(
                int(index) for index in second_generators[second_index]
            )
            projective = exact_backend.recover_projective(
                first_generator,
                second_generator,
                fifth_neighbors,
            )
            if projective is None:
                continue
            denominator_terms = sum(
                bool(coefficient)
                for coefficient in projective.denominator
            )
            numerator_terms = sum(
                bool(coefficient)
                for coefficient in (
                    *projective.numerator_x,
                    *projective.numerator_y,
                )
            )
            coefficient_size = sum(
                abs(int(coefficient.p)).bit_length()
                + int(coefficient.q).bit_length()
                for coefficient in (
                    *projective.numerator_x,
                    *projective.numerator_y,
                    *projective.denominator,
                )
                if coefficient
            )
            certified.append(
                (
                    (
                        denominator_terms,
                        numerator_terms,
                        coefficient_size,
                    ),
                    first_generator,
                    second_generator,
                    projective,
                )
            )
        certified.sort(key=lambda record: record[0])
        for _, first_generator, second_generator, projective in certified:
            candidate = exact_backend.materialize(projective)
            if candidate in existing or candidate in recovered:
                continue
            neighbors = tuple(
                sorted(
                    {
                        *first_generator,
                        *second_generator,
                        projective.fifth_neighbor,
                    }
                )
            )
            neighbor_colors = tuple(
                sorted({int(coloring[index]) for index in neighbors})
            )
            if neighbor_colors != (0, 1, 2, 3, 4):
                continue
            recovered[candidate] = ExactBlockingPoint(
                point=candidate,
                neighbors=neighbors,
                neighbor_colors=neighbor_colors,
                generator_pair=first_generator,
            )
            if max_blockers is not None and len(recovered) >= max_blockers:
                break
        blockers = sorted(
            recovered.values(),
            key=lambda candidate: (
                -len(candidate.neighbors),
                candidate.point.approximate(),
            ),
        )
        return blockers, BlockingSearchStats(
            center_pairs=len(first_center_pairs) + len(second_center_pairs),
            first_partition_pairs=first_pair_count,
            second_partition_pairs=second_pair_count,
            first_distinct_intersections=len(first_keys),
            second_distinct_intersections=len(second_keys),
            common_partition_intersections=common_count,
            nonexisting_common_intersections=len(first_indices),
            fifth_color_numeric_hits=len(numeric_hits),
            exact_blockers=len(blockers),
            elapsed_seconds=time.monotonic() - started,
        )

    if exact_backend is None:
        field, field_points = _field_coordinates(points)
        two = field.convert(2)
        one = field.one
        zero = field.zero
    for hit_index, (first_index, second_index, fifth_neighbors) in enumerate(
        numeric_hits,
        1,
    ):
        first_generator = tuple(
            int(index) for index in first_generators[first_index]
        )
        second_generator = tuple(
            int(index) for index in second_generators[second_index]
        )
        four_neighbors = {
            *first_generator,
            *second_generator,
        }
        ax, ay = field_points[first_generator[0]]
        bx, by = field_points[first_generator[1]]
        cx, cy = field_points[second_generator[0]]
        dx, dy = field_points[second_generator[1]]
        ux, uy = bx - ax, by - ay
        vx, vy = dx - cx, dy - cy
        first_rhs = (
            bx * bx + by * by - ax * ax - ay * ay
        ) / two
        second_rhs = (
            dx * dx + dy * dy - cx * cx - cy * cy
        ) / two
        determinant = ux * vy - uy * vx
        if determinant == zero:
            continue
        candidate_x = (
            first_rhs * vy - uy * second_rhs
        ) / determinant
        candidate_y = (
            ux * second_rhs - first_rhs * vx
        ) / determinant
        if any(
            (candidate_x - field_points[index][0]) ** 2
            + (candidate_y - field_points[index][1]) ** 2
            != one
            for index in four_neighbors
        ):
            continue
        fifth = next(
            (
                index
                for index in fifth_neighbors
                if (
                    (candidate_x - field_points[index][0]) ** 2
                    + (candidate_y - field_points[index][1]) ** 2
                    == one
                )
            ),
            None,
        )
        if fifth is None:
            continue
        candidate = Point(
            field.to_sympy(candidate_x),
            field.to_sympy(candidate_y),
        )
        if candidate in existing:
            continue
        neighbors = tuple(
            sorted(
                {
                    *four_neighbors,
                    fifth,
                }
            )
        )
        neighbor_colors = tuple(
            sorted({int(coloring[index]) for index in neighbors})
        )
        if neighbor_colors != (0, 1, 2, 3, 4):
            continue
        recovered[candidate] = ExactBlockingPoint(
            point=candidate,
            neighbors=neighbors,
            neighbor_colors=neighbor_colors,
            generator_pair=first_generator,
        )
        if max_blockers is not None and len(recovered) >= max_blockers:
            break
        if progress is not None and (
            hit_index % 100 == 0 or hit_index == len(numeric_hits)
        ):
            progress(
                "exact_recovery",
                {
                    "processed": hit_index,
                    "numeric_hits": len(numeric_hits),
                    "exact_blockers": len(recovered),
                },
            )

    blockers = sorted(
        recovered.values(),
        key=lambda candidate: (
            -len(candidate.neighbors),
            candidate.point.approximate(),
        ),
    )
    return blockers, BlockingSearchStats(
        center_pairs=len(first_center_pairs) + len(second_center_pairs),
        first_partition_pairs=first_pair_count,
        second_partition_pairs=second_pair_count,
        first_distinct_intersections=len(first_keys),
        second_distinct_intersections=len(second_keys),
        common_partition_intersections=common_count,
        nonexisting_common_intersections=len(first_indices),
        fifth_color_numeric_hits=len(numeric_hits),
        exact_blockers=len(blockers),
        elapsed_seconds=time.monotonic() - started,
    )
