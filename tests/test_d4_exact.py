import sympy as sp

from d4_exact import (
    D4Algebra,
    D4CoordinateBackend,
    D4ProjectiveCoordinateBackend,
    D4ProjectivePoint,
)
from exact_geometry import Point


def test_d4_algebra_enforces_the_two_quadratic_relations() -> None:
    algebra = D4Algebra(
        sp.Rational(25, 8) - 3 * sp.sqrt(33) / 8,
        sp.Rational(25, 8) + 3 * sp.sqrt(33) / 8,
    )

    assert algebra.mul(algebra.h(0), algebra.h(0)) == algebra.from_base_expr(
        sp.Rational(25, 8) - 3 * sp.sqrt(33) / 8
    )
    assert algebra.mul(algebra.h(1), algebra.h(1)) == algebra.from_base_expr(
        sp.Rational(25, 8) + 3 * sp.sqrt(33) / 8
    )


def test_d4_backend_certifies_a_bisector_intersection_without_division() -> None:
    algebra = D4Algebra(sp.Rational(1, 2), sp.Rational(2))
    points = [
        Point(1, 0),
        Point(-1, 0),
        Point(0, 1),
        Point(0, -1),
        Point(sp.Rational(3, 5), sp.Rational(4, 5)),
    ]
    coordinates = [
        (
            algebra.from_base_expr(point.x),
            algebra.from_base_expr(point.y),
        )
        for point in points
    ]
    backend = D4CoordinateBackend(algebra, coordinates)

    recovered = backend.recover((0, 1), (2, 3), (4,))

    assert recovered is not None
    candidate, fifth = recovered
    assert candidate == Point(0, 0)
    assert fifth == 4
    projective = backend.recover_projective((0, 1), (2, 3), (4,))
    assert projective is not None
    assert backend.unit_neighbors(projective) == (0, 1, 2, 3, 4)
    one = algebra.one
    zero = algebra.zero
    right = D4ProjectivePoint(one, zero, one, -1)
    assert backend.projective_unit(projective, right)
    scaled_origin = D4ProjectivePoint(zero, zero, algebra.rational(2), -1)
    assert backend.projective_equal(projective, scaled_origin)
    four_only = backend.recover_four_projective((0, 1), (2, 3))
    assert four_only is not None
    assert backend.materialize(four_only).approximate() == (0.0, 0.0)


def test_projective_backend_recovers_from_fractional_centers() -> None:
    algebra = D4Algebra(2, 7)

    def center(x_numerator: int, y_numerator: int, denominator: int):
        return D4ProjectivePoint(
            algebra.rational(x_numerator),
            algebra.rational(y_numerator),
            algebra.rational(denominator),
            -1,
        )

    backend = D4ProjectiveCoordinateBackend(
        algebra,
        (
            center(2, 0, 2),
            center(-3, 0, 3),
            center(0, 4, 4),
            center(0, -5, 5),
        ),
    )
    recovered = backend.recover_four_projective((0, 1), (2, 3))

    assert recovered is not None
    assert backend.materialize(recovered) == Point(0, 0)
    assert backend.unit_neighbors(recovered) == (0, 1, 2, 3)


def test_projective_normalization_removes_a_common_rational_scale() -> None:
    algebra = D4Algebra(2, 7)
    scaled = D4ProjectivePoint(
        algebra.scale(algebra.one, sp.Rational(18, 35)),
        algebra.scale(algebra.h(0), sp.Rational(18, 35)),
        algebra.scale(algebra.one, sp.Rational(6, 35)),
        9,
    )

    normalized = algebra.normalize_projective(scaled)

    assert normalized.numerator_x == algebra.rational(3)
    assert normalized.numerator_y == algebra.scale(algebra.h(0), 3)
    assert normalized.denominator == algebra.one
    assert normalized.fifth_neighbor == 9
