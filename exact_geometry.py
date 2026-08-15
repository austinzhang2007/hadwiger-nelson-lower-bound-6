"""Exact algebraic geometry for finite unit-distance graphs.

The certificate-facing API never uses floating-point comparisons.  Numeric
approximations belong in ``cegis.py`` and every retained relation is checked
again here.
"""

from __future__ import annotations

import ast
from functools import lru_cache
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import sympy as sp
from sympy.polys.domains import QQ


RealExpr = sp.Expr
Edge = tuple[int, int]


def _canonical(expr: sp.Expr | int) -> sp.Expr:
    value = sp.sympify(expr)
    return sp.radsimp(sp.cancel(sp.expand(value)))


def _parse_ast(node: ast.AST) -> sp.Expr:
    if isinstance(node, ast.Expression):
        return _parse_ast(node.body)
    if isinstance(node, ast.Constant):
        if type(node.value) is not int:
            raise ValueError("only integer literals are allowed")
        return sp.Integer(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        value = _parse_ast(node.operand)
        return value if isinstance(node.op, ast.UAdd) else -value
    if isinstance(node, ast.BinOp):
        left, right = _parse_ast(node.left), _parse_ast(node.right)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            if right == 0:
                raise ValueError("division by zero")
            return left / right
        raise ValueError("only +, -, *, and / are allowed")
    if (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "sqrt"
        and len(node.args) == 1
        and not node.keywords
    ):
        radicand = _parse_ast(node.args[0])
        if radicand.is_nonnegative is False:
            raise ValueError("coordinates must be real")
        return sp.sqrt(radicand)
    raise ValueError(f"unsupported expression node: {ast.dump(node)}")


def parse_real_expr(text: str) -> sp.Expr:
    """Parse the restricted ``sympy-radical-v1`` expression grammar."""

    if len(text) > 100_000:
        raise ValueError("expression is unreasonably long")
    try:
        tree = ast.parse(text.strip(), mode="eval")
    except SyntaxError as exc:
        raise ValueError(f"invalid algebraic expression: {text!r}") from exc
    value = _canonical(_parse_ast(tree))
    if value.is_real is False:
        raise ValueError("coordinate expression is not real")
    return value


@dataclass(frozen=True, slots=True)
class Point:
    x: sp.Expr
    y: sp.Expr

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", _canonical(self.x))
        object.__setattr__(self, "y", _canonical(self.y))
        if self.x.is_real is False or self.y.is_real is False:
            raise ValueError("Point coordinates must be real")

    @classmethod
    def from_certified_unreduced(
        cls,
        x: sp.Expr,
        y: sp.Expr,
    ) -> "Point":
        """Build a point whose exact field certificate replaces simplification.

        This is reserved for expressions already certified by an exact algebra
        backend.  It avoids potentially explosive radical rationalization while
        retaining the exact symbolic expressions for serialization.
        """

        point = object.__new__(cls)
        object.__setattr__(point, "x", sp.sympify(x))
        object.__setattr__(point, "y", sp.sympify(y))
        return point

    def squared_distance(self, other: "Point") -> sp.Expr:
        dx, dy = self.x - other.x, self.y - other.y
        return _canonical(dx * dx + dy * dy)

    def approximate(self, digits: int = 17) -> tuple[float, float]:
        # Projective D4 coordinates can be quotients of radical sums with
        # severe cancellation.  SymPy's small default ``maxn`` may then stop
        # early and return powers-of-two garbage even though the exact quotient
        # is modest.  A large guard-precision ceiling lets evalf adapt; this
        # remains a search approximation and never certifies an edge.
        return (
            float(sp.N(self.x, digits, maxn=100_000)),
            float(sp.N(self.y, digits, maxn=100_000)),
        )

    def to_json(self, vertex_id: int) -> dict[str, int | str]:
        return {"id": vertex_id, "x": sp.sstr(self.x), "y": sp.sstr(self.y)}


def _split_coordinate_pair(text: str) -> tuple[str, str]:
    depth = 0
    comma = None
    for index, character in enumerate(text):
        if character in "([":
            depth += 1
        elif character in ")]":
            depth -= 1
        elif character == "," and depth == 0:
            if comma is not None:
                raise ValueError(f"too many coordinate components: {text!r}")
            comma = index
    if depth != 0 or comma is None:
        raise ValueError(f"malformed coordinate pair: {text!r}")
    return text[:comma], text[comma + 1 :]


def _mathematica_to_radical(text: str) -> str:
    converted = text.replace("Sqrt[", "sqrt(").replace("]", ")")
    if "[" in converted or "]" in converted or "^" in converted:
        raise ValueError(f"unsupported Mathematica syntax: {text!r}")
    return converted


def parse_mathematica_vertices(path: str | Path) -> list[Point]:
    """Read Heule's author-supplied ``{x, y}`` Mathematica coordinate format."""

    vertices: list[Point] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if not (line.startswith("{") and line.endswith("}")):
                raise ValueError(f"{path}:{line_number}: expected {{x, y}}")
            x_text, y_text = _split_coordinate_pair(line[1:-1])
            vertices.append(
                Point(
                    parse_real_expr(_mathematica_to_radical(x_text)),
                    parse_real_expr(_mathematica_to_radical(y_text)),
                )
            )
    return vertices


def _radical_generators(points: Sequence[Point]) -> list[sp.Expr]:
    radicals: set[sp.Expr] = set()
    for point in points:
        for coordinate in (point.x, point.y):
            for power in coordinate.atoms(sp.Pow):
                if power.exp == sp.Rational(1, 2):
                    radicals.add(power)
    return sorted(radicals, key=sp.default_sort_key)


def _field_coordinates(
    points: Sequence[Point],
) -> tuple[object, list[tuple[object, object]]]:
    generators = _radical_generators(points)
    field = QQ.algebraic_field(*generators) if generators else QQ
    return field, [
        (field.from_sympy(point.x), field.from_sympy(point.y)) for point in points
    ]


@lru_cache(maxsize=16)
def _complete_unit_edges_cached(points: tuple[Point, ...]) -> frozenset[Edge]:
    field, coordinates = _field_coordinates(points)
    one = field.one
    edges: set[Edge] = set()
    for i, (xi, yi) in enumerate(coordinates):
        for j in range(i + 1, len(coordinates)):
            xj, yj = coordinates[j]
            dx, dy = xi - xj, yi - yj
            if dx * dx + dy * dy == one:
                edges.add((i, j))
    return frozenset(edges)


def complete_unit_edges(points: Sequence[Point]) -> set[Edge]:
    """Return every exact unit-distance pair, using zero-based vertex indices."""

    return set(_complete_unit_edges_cached(tuple(points)))


def complete_bipartite_unit_edges(
    left_points: Sequence[Point],
    right_points: Sequence[Point],
) -> set[tuple[int, int]]:
    """Return every exact unit pair with one endpoint in each input sequence."""

    combined = [*left_points, *right_points]
    if not left_points or not right_points:
        return set()
    field, coordinates = _field_coordinates(combined)
    split = len(left_points)
    one = field.one
    edges: set[tuple[int, int]] = set()
    for left_index, (xl, yl) in enumerate(coordinates[:split]):
        for right_index, (xr, yr) in enumerate(coordinates[split:]):
            dx, dy = xl - xr, yl - yr
            if dx * dx + dy * dy == one:
                edges.add((left_index, right_index))
    return edges


def exact_unit_neighbors(point: Point, points: Sequence[Point]) -> tuple[int, ...]:
    combined = [*points, point]
    field, coordinates = _field_coordinates(combined)
    xp, yp = coordinates[-1]
    one = field.one
    neighbors: list[int] = []
    for index, (x, y) in enumerate(coordinates[:-1]):
        dx, dy = xp - x, yp - y
        if dx * dx + dy * dy == one:
            neighbors.append(index)
    return tuple(neighbors)


def _strict_sign(value: sp.Expr) -> int:
    value = _canonical(value)
    if value == 0:
        return 0
    if value.is_positive:
        return 1
    if value.is_negative:
        return -1
    approximation = sp.N(value, 80)
    if approximation > 0:
        return 1
    if approximation < 0:
        return -1
    raise ValueError(f"could not determine exact algebraic sign: {value}")


def unit_circle_intersections(a: Point, b: Point) -> tuple[Point, ...]:
    """Return the exact intersection(s) of two unit circles."""

    dx, dy = b.x - a.x, b.y - a.y
    d2 = _canonical(dx * dx + dy * dy)
    if d2 == 0 or _strict_sign(4 - d2) < 0:
        return ()
    midpoint_x, midpoint_y = (a.x + b.x) / 2, (a.y + b.y) / 2
    if d2 == 4:
        return (Point(midpoint_x, midpoint_y),)
    factor = _canonical(sp.sqrt((4 - d2) / d2) / 2)
    first = Point(midpoint_x - dy * factor, midpoint_y + dx * factor)
    second = Point(midpoint_x + dy * factor, midpoint_y - dx * factor)
    if first == second:
        return (first,)
    return (first, second)


def serialize_graph(
    name: str, points: Sequence[Point], edges: Iterable[Edge]
) -> dict[str, object]:
    return {
        "schema": "udg-exact-v1",
        "name": name,
        "coordinate_format": "sympy-radical-v1",
        "vertices": [
            point.to_json(vertex_id=index + 1) for index, point in enumerate(points)
        ],
        "edges": [[u + 1, v + 1] for u, v in sorted(set(edges))],
    }
