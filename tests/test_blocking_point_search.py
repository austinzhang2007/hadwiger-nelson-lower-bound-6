from exact_geometry import Point
from blocking_point_search import exact_five_color_blockers
from d4_exact import D4Algebra, D4ProjectiveCoordinateBackend, D4ProjectivePoint
import sympy as sp


def test_color_partition_search_recovers_an_exact_five_color_blocker() -> None:
    points = [
        Point(1, 0),
        Point(-1, 0),
        Point(0, 1),
        Point(0, -1),
        Point(sp.Rational(3, 5), sp.Rational(4, 5)),
    ]

    blockers, stats = exact_five_color_blockers(
        points,
        (0, 1, 2, 3, 4),
        quantization_digits=9,
    )

    origin = next(candidate for candidate in blockers if candidate.point == Point(0, 0))
    assert origin.neighbors == (0, 1, 2, 3, 4)
    assert origin.neighbor_colors == (0, 1, 2, 3, 4)
    assert stats.common_partition_intersections >= 1


def test_search_accepts_projective_candidate_centers() -> None:
    algebra = D4Algebra(2, 7)
    points = [
        Point(1, 0),
        Point(-1, 0),
        Point(0, 1),
        Point(0, -1),
        Point(sp.Rational(3, 5), sp.Rational(4, 5)),
    ]

    def projective(point: Point, denominator: int) -> D4ProjectivePoint:
        return D4ProjectivePoint(
            algebra.from_base_expr(point.x * denominator),
            algebra.from_base_expr(point.y * denominator),
            algebra.rational(denominator),
            -1,
        )

    backend = D4ProjectiveCoordinateBackend(
        algebra,
        tuple(
            projective(point, denominator)
            for point, denominator in zip(points, (2, 3, 4, 5, 6))
        ),
    )
    blockers, stats = exact_five_color_blockers(
        points,
        (0, 1, 2, 3, 4),
        exact_backend=backend,
    )

    assert any(blocker.point == Point(0, 0) for blocker in blockers)
    assert stats.exact_blockers >= 1
