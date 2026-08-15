from pathlib import Path

import pytest

from coloring_sat import ColoringSAT, parse_dimacs_edge
from tests.test_exact_geometry import moser_spindle_points
from exact_geometry import complete_unit_edges


ROOT = Path(__file__).resolve().parents[1]


def test_moser_spindle_is_not_three_colorable_but_is_four_colorable() -> None:
    edges = complete_unit_edges(moser_spindle_points())
    assert ColoringSAT(7, edges, 3).solve() is None
    coloring = ColoringSAT(7, edges, 4).solve()
    assert coloring is not None
    assert len(coloring) == 7


def test_color_precedence_enumeration_removes_color_permutations() -> None:
    solver = ColoringSAT(3, {(0, 1), (1, 2), (0, 2)}, 3)
    models = list(solver.enumerate_colorings(limit=10))
    assert models == [(0, 1, 2)]


def test_public_heule_graph_is_four_unsat_and_five_sat() -> None:
    n, edges = parse_dimacs_edge(ROOT / "third_party/CNP-SAT/edge/529.edge")
    assert n == 529
    assert len(edges) == 2670
    assert ColoringSAT(n, edges, 4).solve() is None
    assert ColoringSAT(n, edges, 5).solve() is not None


def test_rejects_edges_outside_the_graph() -> None:
    with pytest.raises(ValueError, match="outside graph"):
        ColoringSAT(2, {(0, 2)}, 3)


def test_symmetry_clique_fixes_distinct_colors_with_unit_clauses() -> None:
    triangle = {(0, 1), (1, 2), (0, 2)}
    encoding = ColoringSAT(
        4,
        triangle,
        4,
        break_color_symmetry=False,
        symmetry_clique=(0, 1, 2),
    )

    coloring = encoding.solve()

    assert coloring is not None
    assert coloring[:3] == (0, 1, 2)
    assert all(
        [encoding.variable(vertex, vertex)] in encoding.clauses
        for vertex in range(3)
    )


def test_symmetry_clique_rejects_a_nonedge() -> None:
    with pytest.raises(ValueError, match="must induce a clique"):
        ColoringSAT(
            3,
            {(0, 1)},
            3,
            break_color_symmetry=False,
            symmetry_clique=(0, 1, 2),
        )


def test_solve_accepts_a_complete_phase_hint() -> None:
    coloring = ColoringSAT(
        2,
        {(0, 1)},
        2,
        break_color_symmetry=False,
    ).solve(phases=(1, -2, -3, 4))

    assert coloring == (0, 1)
