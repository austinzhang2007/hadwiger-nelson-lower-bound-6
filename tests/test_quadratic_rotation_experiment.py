import sympy as sp

from exact_geometry import Point


def test_new_rotation_is_unit_and_contains_sqrt7():
    from quadratic_rotation_experiment import rotate_sqrt7
    p = rotate_sqrt7(Point(1, 0))
    assert p.x == sp.Rational(1, 8)
    assert p.y == 3 * sp.sqrt(7) / 8
    assert sp.expand(p.x**2 + p.y**2) == 1


def test_integer_complete_edges_agree_with_exact_geometry():
    from quadratic_rotation_experiment import integer_complete_edges
    from exact_geometry import complete_unit_edges
    points = [Point(0, 0), Point(1, 0), Point(sp.Rational(1, 2), sp.sqrt(3)/2),
              Point(sp.Rational(1, 8), 3*sp.sqrt(7)/8),
              Point(1 + sp.Rational(1, 10**20), 0)]
    assert integer_complete_edges(points) == complete_unit_edges(points)


def test_integer_complete_edges_rejects_unrepresented_radical():
    import pytest
    from quadratic_rotation_experiment import integer_complete_edges
    with pytest.raises(ValueError, match="field"):
        integer_complete_edges([Point(0, 0), Point(sp.sqrt(13), 0)])
