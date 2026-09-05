from pysat.solvers import Solver

from coloring_sat import ColoringSAT
from lean_coloring_cnf import (
    decode_lean_coloring_model,
    lean_coloring_clauses,
)


def solve_lean(vertex_count, edges, colors):
    clauses = lean_coloring_clauses(vertex_count, edges, colors)
    with Solver(name="glucose4", bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None
        return decode_lean_coloring_model(
            solver.get_model(),
            vertex_count=vertex_count,
            edges=edges,
            colors=colors,
        )


def test_lean_encoding_decodes_multicolor_model_to_proper_coloring():
    coloring = decode_lean_coloring_model(
        [1, 2, -3, -4, -5, 6],
        vertex_count=2,
        edges={(0, 1)},
        colors=3,
    )

    assert coloring == (0, 2)


def test_lean_encoding_matches_exactly_one_on_sat_and_unsat_graphs():
    instances = [
        (3, {(0, 1), (1, 2)}, 2),
        (3, {(0, 1), (1, 2), (0, 2)}, 2),
        (4, {(0, 1), (1, 2), (2, 3), (3, 0)}, 2),
    ]

    for vertex_count, edges, colors in instances:
        lean = solve_lean(vertex_count, edges, colors)
        exact = ColoringSAT(
            vertex_count,
            edges,
            colors,
            break_color_symmetry=False,
        ).solve()
        assert (lean is None) == (exact is None)


def test_symmetry_clique_is_fixed_without_at_most_one_clauses():
    clauses = lean_coloring_clauses(
        3,
        {(0, 1), (1, 2), (0, 2)},
        3,
        symmetry_clique=(0, 1, 2),
    )

    assert [1] in clauses
    assert [5] in clauses
    assert [9] in clauses
    assert len(clauses) == 3 + 3 * 3 + 3
