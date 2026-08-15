from pysat.formula import CNF
from pysat.solvers import Solver

from coloring_sat import ColoringSAT
from residual_symmetry_cnf import add_residual_color_precedence


def test_residual_precedence_uses_linear_auxiliaries(tmp_path) -> None:
    source = tmp_path / "one.cnf"
    target = tmp_path / "one-sym.cnf"
    ColoringSAT(
        1,
        (),
        5,
        break_color_symmetry=False,
    ).write_dimacs_cnf(source)

    stats = add_residual_color_precedence(source, target, vertex_count=1)
    formula = CNF(from_file=str(target))

    assert stats["added_variables"] == 1
    assert stats["added_clauses"] == 3
    with Solver(name="cadical195", bootstrap_with=formula.clauses) as solver:
        assert solver.solve(assumptions=[4])
        assert not solver.solve(assumptions=[5])
