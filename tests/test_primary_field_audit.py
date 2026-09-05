import random

import pytest
import sympy as sp

from d4_exact import D4Algebra


def test_independent_multiplication_matches_all_basis_products():
    from primary_field_audit import PrimaryField

    field = PrimaryField()
    reference = D4Algebra((25 - 3 * sp.sqrt(33)) / 8,
                          (25 + 3 * sp.sqrt(33)) / 8)
    for i in range(32):
        for j in range(32):
            a = tuple(int(k == i) for k in range(32))
            b = tuple(int(k == j) for k in range(32))
            assert tuple(map(str, field.mul(a, b))) == tuple(
                map(str, reference.mul(a, b))
            )


def test_dense_rational_products_match_sympy_reference():
    from primary_field_audit import PrimaryField

    rng = random.Random(20260906)
    field = PrimaryField()
    reference = D4Algebra((25 - 3 * sp.sqrt(33)) / 8,
                          (25 + 3 * sp.sqrt(33)) / 8)
    for _ in range(4):
        a = tuple(sp.Rational(rng.randrange(-50, 51), rng.randrange(1, 20))
                  for _ in range(32))
        b = tuple(sp.Rational(rng.randrange(-50, 51), rng.randrange(1, 20))
                  for _ in range(32))
        assert tuple(map(str, field.mul(field.element(a), field.element(b)))) == tuple(
            map(str, reference.mul(a, b))
        )


def test_unit_identity_rejects_zero_denominator_and_near_unit_distance():
    from primary_field_audit import PrimaryField

    f = PrimaryField()
    origin = (f.zero, f.zero, f.one)
    unit = (f.one, f.zero, f.one)
    assert f.unit(origin, unit)
    almost = (f.element(["100000000000000000001/100000000000000000000"] + [0] * 31),
              f.zero, f.one)
    assert not f.unit(origin, almost)
    with pytest.raises(ValueError, match="denominator"):
        f.unit(origin, (f.zero, f.zero, f.zero))


def test_homogeneous_unit_identity_is_scale_invariant():
    from primary_field_audit import PrimaryField

    f = PrimaryField()
    r = lambda n: f.element([n] + [0] * 31)
    assert f.unit((r(0), r(0), r(7)), (r(3), r(4), r(5)))
    assert not f.unit((r(0), r(0), r(7)), (r(3), r(4), r(6)))
