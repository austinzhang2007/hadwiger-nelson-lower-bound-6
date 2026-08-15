from pathlib import Path

import sympy as sp

from exact_geometry import (
    Point,
    complete_bipartite_unit_edges,
    complete_unit_edges,
    parse_mathematica_vertices,
    parse_real_expr,
    unit_circle_intersections,
)


ROOT = Path(__file__).resolve().parents[1]


def moser_spindle_points() -> list[Point]:
    s3, s11, s33 = sp.sqrt(3), sp.sqrt(11), sp.sqrt(33)
    return [
        Point(sp.Integer(0), sp.Integer(0)),
        Point(s3, sp.Integer(0)),
        Point(s3 / 2, sp.Rational(1, 2)),
        Point(s3 / 2, sp.Rational(-1, 2)),
        Point(5 * s3 / 6, s33 / 6),
        Point((5 * s3 - s11) / 12, (s33 + 5) / 12),
        Point((5 * s3 + s11) / 12, (s33 - 5) / 12),
    ]


def test_restricted_expression_parser_is_exact() -> None:
    value = parse_real_expr("(3 + sqrt(33))/6")
    assert sp.expand(value * 6 - 3) ** 2 == 33


def test_certified_unreduced_point_preserves_an_exact_expression() -> None:
    expression = sp.Mul(
        sp.sqrt(3) + 1,
        sp.Pow(sp.sqrt(3) + 1, -1, evaluate=False),
        evaluate=False,
    )

    point = Point.from_certified_unreduced(expression, sp.S.Zero)

    assert point.x is expression
    assert point.approximate() == (1.0, 0.0)


def test_moser_spindle_has_exactly_eleven_unit_edges() -> None:
    edges = complete_unit_edges(moser_spindle_points())
    assert len(edges) == 11
    assert (1, 4) in edges


def test_complete_bipartite_unit_edges_scans_every_cross_pair() -> None:
    left = [Point(0, 0), Point(2, 0)]
    right = [Point(1, 0), Point(0, 1)]

    assert complete_bipartite_unit_edges(left, right) == {
        (0, 0),
        (1, 0),
        (0, 1),
    }


def test_unit_circle_intersections_are_symbolically_unit_distance() -> None:
    a = Point(sp.Integer(0), sp.Integer(0))
    b = Point(sp.Integer(1), sp.Integer(0))
    intersections = unit_circle_intersections(a, b)
    assert len(intersections) == 2
    assert {sp.simplify(p.y) for p in intersections} == {
        -sp.sqrt(3) / 2,
        sp.sqrt(3) / 2,
    }
    for point in intersections:
        assert point.squared_distance(a) == 1
        assert point.squared_distance(b) == 1


def test_author_heule_529_coordinates_and_complete_edges() -> None:
    vertices = parse_mathematica_vertices(
        ROOT / "third_party/CNP-SAT/vtx/529.vtx"
    )
    assert len(vertices) == 529
    assert len(set(vertices)) == 529
    assert len(complete_unit_edges(vertices)) == 2670
