import sympy as sp

from d4_exact import D4Algebra, D4CoordinateBackend
from d4_pair_search import exact_four_color_forced_points
from exact_geometry import Point


def test_four_color_partition_recovers_a_forced_point() -> None:
    points = [
        Point(1, 0),
        Point(-1, 0),
        Point(0, 1),
        Point(0, -1),
        Point(10, 10),
    ]
    algebra = D4Algebra(sp.Rational(1, 2), sp.Rational(2))
    coordinates = [
        (
            algebra.from_base_expr(point.x),
            algebra.from_base_expr(point.y),
        )
        for point in points
    ]

    candidates, stats = exact_four_color_forced_points(
        points,
        (0, 1, 2, 3, 4),
        D4CoordinateBackend(algebra, coordinates),
        missing_color=4,
    )

    origin = next(item for item in candidates if item.point.approximate() == (0.0, 0.0))
    assert origin.neighbors == (0, 1, 2, 3)
    assert origin.forced_color == 4
    assert stats.exact_candidates >= 1
