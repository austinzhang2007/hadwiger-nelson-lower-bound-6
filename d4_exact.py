"""Small exact algebra for the two primary-D4 translation extensions.

The basis is
``sqrt(3)^a sqrt(5)^b sqrt(11)^c h_0^d h_1^e``
with binary exponents.  Arithmetic therefore stays in 32 rational
coefficients and avoids repeated primitive-element factorization in SymPy.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, lcm
import sys
from typing import Sequence

import sympy as sp

from exact_geometry import Point


if hasattr(sys, "set_int_max_str_digits"):
    # Exact projective certificates can legitimately exceed Python's 4300
    # decimal-digit default before homogeneous normalization is applied.
    sys.set_int_max_str_digits(1_000_000)


Element = tuple[sp.Rational, ...]


def _zero_coefficients() -> list[sp.Rational]:
    return [sp.S.Zero] * 32


class D4Algebra:
    """Exact 32-dimensional algebra containing both quadratic translations."""

    _base_squares = ((1, 3), (2, 5), (4, 11))

    def __init__(self, first_h_squared: sp.Expr, second_h_squared: sp.Expr):
        self.zero: Element = tuple(_zero_coefficients())
        self.one = self.rational(1)
        self._radicals = (
            sp.sqrt(3),
            sp.sqrt(5),
            sp.sqrt(11),
            sp.sqrt(first_h_squared),
            sp.sqrt(second_h_squared),
        )
        self._h_square_relations = (
            self.from_base_expr(first_h_squared),
            self.from_base_expr(second_h_squared),
        )

    def rational(self, value: int | sp.Rational) -> Element:
        coefficients = _zero_coefficients()
        coefficients[0] = sp.Rational(value)
        return tuple(coefficients)

    def h(self, kind: int) -> Element:
        if kind not in (0, 1):
            raise ValueError("D4 extension kind must be 0 or 1")
        coefficients = _zero_coefficients()
        coefficients[1 << (3 + kind)] = sp.S.One
        return tuple(coefficients)

    def from_base_expr(self, expression: sp.Expr | int) -> Element:
        s3, s5, s11 = sp.symbols("s3 s5 s11")
        replacements = {
            sp.sqrt(165): s3 * s5 * s11,
            sp.sqrt(55): s5 * s11,
            sp.sqrt(33): s3 * s11,
            sp.sqrt(15): s3 * s5,
            sp.sqrt(11): s11,
            sp.sqrt(5): s5,
            sp.sqrt(3): s3,
        }
        replaced = sp.expand(sp.sympify(expression).xreplace(replacements))
        polynomial = sp.Poly(replaced, s3, s5, s11, domain=sp.QQ)
        coefficients = _zero_coefficients()
        for monomial, coefficient in polynomial.terms():
            if any(exponent not in (0, 1) for exponent in monomial):
                raise ValueError("expression is outside Q(sqrt(3),sqrt(5),sqrt(11))")
            mask = sum(
                (1 << bit) for bit, exponent in enumerate(monomial) if exponent
            )
            coefficients[mask] += sp.Rational(coefficient)
        return tuple(coefficients)

    @staticmethod
    def add(left: Element, right: Element) -> Element:
        return tuple(a + b for a, b in zip(left, right))

    @staticmethod
    def sub(left: Element, right: Element) -> Element:
        return tuple(a - b for a, b in zip(left, right))

    @staticmethod
    def scale(value: Element, scalar: int | sp.Rational) -> Element:
        factor = sp.Rational(scalar)
        return tuple(factor * coefficient for coefficient in value)

    def _multiply_by_base_mask(
        self,
        terms: dict[int, sp.Rational],
        base_mask: int,
        relation_coefficient: sp.Rational,
    ) -> dict[int, sp.Rational]:
        result: dict[int, sp.Rational] = {}
        for mask, coefficient in terms.items():
            common = mask & base_mask
            factor = relation_coefficient
            for bit, square in self._base_squares:
                if common & bit:
                    factor *= square
            output_mask = mask ^ base_mask
            result[output_mask] = (
                result.get(output_mask, sp.S.Zero) + coefficient * factor
            )
        return result

    def mul(self, left: Element, right: Element) -> Element:
        output = _zero_coefficients()
        left_terms = [(i, c) for i, c in enumerate(left) if c]
        right_terms = [(i, c) for i, c in enumerate(right) if c]
        for left_mask, left_coefficient in left_terms:
            for right_mask, right_coefficient in right_terms:
                common = left_mask & right_mask
                scalar = left_coefficient * right_coefficient
                for bit, square in self._base_squares:
                    if common & bit:
                        scalar *= square
                terms = {left_mask ^ right_mask: scalar}
                for kind, h_bit in enumerate((8, 16)):
                    if not common & h_bit:
                        continue
                    expanded: dict[int, sp.Rational] = {}
                    relation = self._h_square_relations[kind]
                    for relation_mask, relation_coefficient in enumerate(relation):
                        if not relation_coefficient:
                            continue
                        contribution = self._multiply_by_base_mask(
                            terms,
                            relation_mask,
                            relation_coefficient,
                        )
                        for mask, coefficient in contribution.items():
                            expanded[mask] = (
                                expanded.get(mask, sp.S.Zero) + coefficient
                            )
                    terms = expanded
                for mask, coefficient in terms.items():
                    output[mask] += coefficient
        return tuple(output)

    def square(self, value: Element) -> Element:
        return self.mul(value, value)

    def to_sympy(self, value: Element) -> sp.Expr:
        return sp.factor(self.to_sympy_unfactored(value))

    def to_sympy_unfactored(self, value: Element) -> sp.Expr:
        """Translate coefficients without invoking polynomial factorization."""

        expression = sp.S.Zero
        for mask, coefficient in enumerate(value):
            if not coefficient:
                continue
            term: sp.Expr = coefficient
            for bit, radical in enumerate(self._radicals):
                if mask & (1 << bit):
                    term *= radical
            expression += term
        return expression

    def projective_coordinate_expr(
        self,
        numerator: Element,
        denominator: Element,
    ) -> sp.Expr:
        """Build a certified quotient without expensive radical factoring."""

        if numerator == self.zero:
            return sp.S.Zero
        raw_numerator = self.to_sympy_unfactored(numerator)
        if denominator == self.one:
            return raw_numerator
        raw_denominator = self.to_sympy_unfactored(denominator)
        return sp.Mul(
            raw_numerator,
            sp.Pow(raw_denominator, -1, evaluate=False),
            evaluate=False,
        )

    def normalize_projective(
        self,
        point: "D4ProjectivePoint",
    ) -> "D4ProjectivePoint":
        """Remove the common rational scale from 96 homogeneous coefficients."""

        values = (
            *point.numerator_x,
            *point.numerator_y,
            *point.denominator,
        )
        common_denominator = 1
        for value in values:
            if value:
                common_denominator = lcm(
                    common_denominator,
                    int(value.q),
                )
        integers = [
            int(value.p) * (common_denominator // int(value.q))
            for value in values
        ]
        common_factor = 0
        for value in integers:
            common_factor = gcd(common_factor, abs(value))
        if common_factor == 0:
            raise ValueError("zero projective coordinate")
        integers = [value // common_factor for value in integers]
        first_nonzero = next(
            value
            for value in (
                *integers[64:96],
                *integers[:64],
            )
            if value
        )
        if first_nonzero < 0:
            integers = [-value for value in integers]
        return D4ProjectivePoint(
            tuple(sp.Rational(value) for value in integers[:32]),
            tuple(sp.Rational(value) for value in integers[32:64]),
            tuple(sp.Rational(value) for value in integers[64:96]),
            point.fifth_neighbor,
        )


@dataclass(frozen=True, slots=True)
class D4ProjectivePoint:
    numerator_x: Element
    numerator_y: Element
    denominator: Element
    fifth_neighbor: int


@dataclass(frozen=True, slots=True)
class D4CoordinateBackend:
    algebra: D4Algebra
    coordinates: Sequence[tuple[Element, Element]]

    def _unit_numerator(
        self,
        numerator_x: Element,
        numerator_y: Element,
        denominator: Element,
        point_index: int,
    ) -> bool:
        x, y = self.coordinates[point_index]
        dx = self.algebra.sub(
            numerator_x,
            self.algebra.mul(x, denominator),
        )
        dy = self.algebra.sub(
            numerator_y,
            self.algebra.mul(y, denominator),
        )
        return self.algebra.add(
            self.algebra.square(dx),
            self.algebra.square(dy),
        ) == self.algebra.square(denominator)

    def recover(
        self,
        first_generator: tuple[int, int],
        second_generator: tuple[int, int],
        fifth_neighbors: Sequence[int],
    ) -> tuple[Point, int] | None:
        """Certify a four-bisector candidate and one fifth-color unit edge."""

        recovered = self.recover_projective(
            first_generator,
            second_generator,
            fifth_neighbors,
        )
        if recovered is None:
            return None
        return self.materialize(recovered), recovered.fifth_neighbor

    def recover_projective(
        self,
        first_generator: tuple[int, int],
        second_generator: tuple[int, int],
        fifth_neighbors: Sequence[int],
    ) -> D4ProjectivePoint | None:
        """Recover a candidate as ``(Nx/D, Ny/D)`` without field division."""

        recovered = self.recover_four_projective(
            first_generator,
            second_generator,
        )
        if recovered is None:
            return None
        fifth = next(
            (
                index
                for index in fifth_neighbors
                if self._unit_numerator(
                    recovered.numerator_x,
                    recovered.numerator_y,
                    recovered.denominator,
                    index,
                )
            ),
            None,
        )
        if fifth is None:
            return None
        return D4ProjectivePoint(
            numerator_x=recovered.numerator_x,
            numerator_y=recovered.numerator_y,
            denominator=recovered.denominator,
            fifth_neighbor=fifth,
        )

    def recover_four_projective(
        self,
        first_generator: tuple[int, int],
        second_generator: tuple[int, int],
    ) -> D4ProjectivePoint | None:
        """Recover and certify a point unit from four generator centers."""

        ax, ay = self.coordinates[first_generator[0]]
        bx, by = self.coordinates[first_generator[1]]
        cx, cy = self.coordinates[second_generator[0]]
        dx, dy = self.coordinates[second_generator[1]]
        ux, uy = self.algebra.sub(bx, ax), self.algebra.sub(by, ay)
        vx, vy = self.algebra.sub(dx, cx), self.algebra.sub(dy, cy)
        first_rhs = self.algebra.scale(
            self.algebra.sub(
                self.algebra.add(
                    self.algebra.square(bx),
                    self.algebra.square(by),
                ),
                self.algebra.add(
                    self.algebra.square(ax),
                    self.algebra.square(ay),
                ),
            ),
            sp.Rational(1, 2),
        )
        second_rhs = self.algebra.scale(
            self.algebra.sub(
                self.algebra.add(
                    self.algebra.square(dx),
                    self.algebra.square(dy),
                ),
                self.algebra.add(
                    self.algebra.square(cx),
                    self.algebra.square(cy),
                ),
            ),
            sp.Rational(1, 2),
        )
        determinant = self.algebra.sub(
            self.algebra.mul(ux, vy),
            self.algebra.mul(uy, vx),
        )
        if determinant == self.algebra.zero:
            return None
        numerator_x = self.algebra.sub(
            self.algebra.mul(first_rhs, vy),
            self.algebra.mul(uy, second_rhs),
        )
        numerator_y = self.algebra.sub(
            self.algebra.mul(ux, second_rhs),
            self.algebra.mul(first_rhs, vx),
        )
        four_neighbors = {
            *first_generator,
            *second_generator,
        }
        recovered = self.algebra.normalize_projective(
            D4ProjectivePoint(
                numerator_x=numerator_x,
                numerator_y=numerator_y,
                denominator=determinant,
                fifth_neighbor=-1,
            )
        )
        if any(
            not self._unit_numerator(
                recovered.numerator_x,
                recovered.numerator_y,
                recovered.denominator,
                index,
            )
            for index in four_neighbors
        ):
            return None
        return recovered

    def materialize(self, point: D4ProjectivePoint) -> Point:
        return Point.from_certified_unreduced(
            self.algebra.projective_coordinate_expr(
                point.numerator_x,
                point.denominator,
            ),
            self.algebra.projective_coordinate_expr(
                point.numerator_y,
                point.denominator,
            ),
        )

    def unit_neighbors(self, point: D4ProjectivePoint) -> tuple[int, ...]:
        """Enumerate every unit neighbor in the reconstructed exact graph."""

        return tuple(
            index
            for index in range(len(self.coordinates))
            if self._unit_numerator(
                point.numerator_x,
                point.numerator_y,
                point.denominator,
                index,
            )
        )

    def projective_unit(
        self,
        left: D4ProjectivePoint,
        right: D4ProjectivePoint,
    ) -> bool:
        """Check the unit relation between two projective algebra points."""

        dx = self.algebra.sub(
            self.algebra.mul(left.numerator_x, right.denominator),
            self.algebra.mul(right.numerator_x, left.denominator),
        )
        dy = self.algebra.sub(
            self.algebra.mul(left.numerator_y, right.denominator),
            self.algebra.mul(right.numerator_y, left.denominator),
        )
        denominator = self.algebra.mul(
            left.denominator,
            right.denominator,
        )
        return self.algebra.add(
            self.algebra.square(dx),
            self.algebra.square(dy),
        ) == self.algebra.square(denominator)

    def projective_equal(
        self,
        left: D4ProjectivePoint,
        right: D4ProjectivePoint,
    ) -> bool:
        """Check equality without choosing a normalized denominator."""

        return self.algebra.mul(
            left.numerator_x,
            right.denominator,
        ) == self.algebra.mul(
            right.numerator_x,
            left.denominator,
        ) and self.algebra.mul(
            left.numerator_y,
            right.denominator,
        ) == self.algebra.mul(
            right.numerator_y,
            left.denominator,
        )

    def as_projective_backend(
        self,
        extra: Sequence[D4ProjectivePoint] = (),
    ) -> "D4ProjectiveCoordinateBackend":
        """Promote affine field coordinates and append rational field points."""

        one = self.algebra.one
        return D4ProjectiveCoordinateBackend(
            self.algebra,
            tuple(
                D4ProjectivePoint(x, y, one, -1)
                for x, y in self.coordinates
            )
            + tuple(self.algebra.normalize_projective(point) for point in extra),
        )


@dataclass(frozen=True, slots=True)
class D4ProjectiveCoordinateBackend:
    """Exact backend whose circle centers may themselves have denominators.

    A bisector equation is cleared of all center denominators before Cramer's
    rule is applied.  Returned coordinates therefore remain projective
    elements of the same 32-dimensional algebra; no floating-point decision
    enters recovery or unit-distance certification.
    """

    algebra: D4Algebra
    coordinates: Sequence[D4ProjectivePoint]

    def _unit_numerator(
        self,
        numerator_x: Element,
        numerator_y: Element,
        denominator: Element,
        point_index: int,
    ) -> bool:
        point = self.coordinates[point_index]
        dx = self.algebra.sub(
            self.algebra.mul(numerator_x, point.denominator),
            self.algebra.mul(point.numerator_x, denominator),
        )
        dy = self.algebra.sub(
            self.algebra.mul(numerator_y, point.denominator),
            self.algebra.mul(point.numerator_y, denominator),
        )
        common_denominator = self.algebra.mul(
            denominator,
            point.denominator,
        )
        return self.algebra.add(
            self.algebra.square(dx),
            self.algebra.square(dy),
        ) == self.algebra.square(common_denominator)

    def _bisector_row(
        self,
        left_index: int,
        right_index: int,
    ) -> tuple[Element, Element, Element]:
        left = self.coordinates[left_index]
        right = self.coordinates[right_index]
        common_denominator = self.algebra.mul(
            left.denominator,
            right.denominator,
        )
        ux = self.algebra.sub(
            self.algebra.mul(right.numerator_x, left.denominator),
            self.algebra.mul(left.numerator_x, right.denominator),
        )
        uy = self.algebra.sub(
            self.algebra.mul(right.numerator_y, left.denominator),
            self.algebra.mul(left.numerator_y, right.denominator),
        )
        left_norm = self.algebra.add(
            self.algebra.square(left.numerator_x),
            self.algebra.square(left.numerator_y),
        )
        right_norm = self.algebra.add(
            self.algebra.square(right.numerator_x),
            self.algebra.square(right.numerator_y),
        )
        rhs = self.algebra.sub(
            self.algebra.mul(
                right_norm,
                self.algebra.square(left.denominator),
            ),
            self.algebra.mul(
                left_norm,
                self.algebra.square(right.denominator),
            ),
        )
        twice_common = self.algebra.scale(common_denominator, 2)
        return (
            self.algebra.mul(twice_common, ux),
            self.algebra.mul(twice_common, uy),
            rhs,
        )

    def recover_four_projective(
        self,
        first_generator: tuple[int, int],
        second_generator: tuple[int, int],
    ) -> D4ProjectivePoint | None:
        """Recover a four-center candidate with denominator clearing."""

        ax, ay, arhs = self._bisector_row(*first_generator)
        bx, by, brhs = self._bisector_row(*second_generator)
        determinant = self.algebra.sub(
            self.algebra.mul(ax, by),
            self.algebra.mul(ay, bx),
        )
        if determinant == self.algebra.zero:
            return None
        numerator_x = self.algebra.sub(
            self.algebra.mul(arhs, by),
            self.algebra.mul(ay, brhs),
        )
        numerator_y = self.algebra.sub(
            self.algebra.mul(ax, brhs),
            self.algebra.mul(arhs, bx),
        )
        recovered = self.algebra.normalize_projective(
            D4ProjectivePoint(
                numerator_x,
                numerator_y,
                determinant,
                -1,
            )
        )
        if any(
            not self._unit_numerator(
                recovered.numerator_x,
                recovered.numerator_y,
                recovered.denominator,
                index,
            )
            for index in {*first_generator, *second_generator}
        ):
            return None
        return recovered

    def recover_projective(
        self,
        first_generator: tuple[int, int],
        second_generator: tuple[int, int],
        fifth_neighbors: Sequence[int],
    ) -> D4ProjectivePoint | None:
        recovered = self.recover_four_projective(
            first_generator,
            second_generator,
        )
        if recovered is None:
            return None
        fifth = next(
            (
                index
                for index in fifth_neighbors
                if self._unit_numerator(
                    recovered.numerator_x,
                    recovered.numerator_y,
                    recovered.denominator,
                    index,
                )
            ),
            None,
        )
        if fifth is None:
            return None
        return D4ProjectivePoint(
            recovered.numerator_x,
            recovered.numerator_y,
            recovered.denominator,
            fifth,
        )

    def materialize(self, point: D4ProjectivePoint) -> Point:
        return Point.from_certified_unreduced(
            self.algebra.projective_coordinate_expr(
                point.numerator_x,
                point.denominator,
            ),
            self.algebra.projective_coordinate_expr(
                point.numerator_y,
                point.denominator,
            ),
        )

    def unit_neighbors(self, point: D4ProjectivePoint) -> tuple[int, ...]:
        return tuple(
            index
            for index in range(len(self.coordinates))
            if self._unit_numerator(
                point.numerator_x,
                point.numerator_y,
                point.denominator,
                index,
            )
        )

    def projective_unit(
        self,
        left: D4ProjectivePoint,
        right: D4ProjectivePoint,
    ) -> bool:
        dx = self.algebra.sub(
            self.algebra.mul(left.numerator_x, right.denominator),
            self.algebra.mul(right.numerator_x, left.denominator),
        )
        dy = self.algebra.sub(
            self.algebra.mul(left.numerator_y, right.denominator),
            self.algebra.mul(right.numerator_y, left.denominator),
        )
        denominator = self.algebra.mul(
            left.denominator,
            right.denominator,
        )
        return self.algebra.add(
            self.algebra.square(dx),
            self.algebra.square(dy),
        ) == self.algebra.square(denominator)

    def projective_equal(
        self,
        left: D4ProjectivePoint,
        right: D4ProjectivePoint,
    ) -> bool:
        return self.algebra.mul(
            left.numerator_x,
            right.denominator,
        ) == self.algebra.mul(
            right.numerator_x,
            left.denominator,
        ) and self.algebra.mul(
            left.numerator_y,
            right.denominator,
        ) == self.algebra.mul(
            right.numerator_y,
            left.denominator,
        )


def primary_d4_coordinate_backend(
    base_points: Sequence[Point],
    source_points: Sequence[Point],
    selected_r_indices: Sequence[int],
) -> D4CoordinateBackend:
    """Reconstruct the exact 2510-plus-copy coordinates in the 32-term basis."""

    def canonical(value: sp.Expr) -> sp.Expr:
        return sp.radsimp(sp.cancel(sp.expand(value)))

    h_squares: list[sp.Expr] = []
    selected_h_squares: list[sp.Expr] = []
    for source_index in selected_r_indices:
        r = source_points[source_index]
        squared = canonical(r.x * r.x + r.y * r.y)
        h_squared = canonical((4 - squared) / (4 * squared))
        selected_h_squares.append(h_squared)
        if not any(canonical(h_squared - known) == 0 for known in h_squares):
            h_squares.append(h_squared)
    if len(h_squares) != 2:
        raise ValueError("primary D4 construction must use two conjugate h-squares")
    h_squares.sort(key=lambda value: float(sp.N(value)))
    algebra = D4Algebra(h_squares[0], h_squares[1])

    def base_coordinate(point: Point) -> tuple[Element, Element]:
        return (
            algebra.from_base_expr(point.x),
            algebra.from_base_expr(point.y),
        )

    coordinates = [base_coordinate(point) for point in base_points]
    source_coordinates = [base_coordinate(point) for point in source_points]
    for selected_index, h_squared in zip(
        selected_r_indices,
        selected_h_squares,
    ):
        kind = next(
            index
            for index, known in enumerate(h_squares)
            if canonical(h_squared - known) == 0
        )
        rx, ry = source_coordinates[selected_index]
        h = algebra.h(kind)
        translation_x = algebra.add(
            algebra.scale(rx, sp.Rational(1, 2)),
            algebra.mul(ry, h),
        )
        translation_y = algebra.sub(
            algebra.scale(ry, sp.Rational(1, 2)),
            algebra.mul(rx, h),
        )
        coordinates.extend(
            (
                algebra.add(x, translation_x),
                algebra.add(y, translation_y),
            )
            for x, y in source_coordinates
        )
    return D4CoordinateBackend(algebra, coordinates)
