import sympy as sp

from exact_geometry import Point


def test_anchor_angles_force_unit_chords_in_new_prime_fields():
    from quadratic_rotation_experiment import rotate_exact, integer_complete_edges
    for q, c, s, primes in (
        (sp.Rational(5, 3), sp.Rational(7, 10), sp.sqrt(51)/10, (3, 5, 11, 17)),
        (sp.Rational(4, 3), sp.Rational(5, 8), sp.sqrt(39)/8, (3, 5, 11, 13)),
    ):
        p = Point(sp.sqrt(q), 0)
        r = rotate_exact(p, c, s)
        assert sp.simplify((p.x-r.x)**2 + (p.y-r.y)**2) == 1
        assert integer_complete_edges([p, r], primes=primes) == {(0, 1)}
