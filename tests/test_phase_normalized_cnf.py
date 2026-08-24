from pathlib import Path

from phase_normalized_cnf import (
    canonicalize_residual_coloring,
    phase_true_variables,
    rewrite_cnf_for_false_phase,
    unflip_model,
)


def test_canonicalize_swaps_residual_colors_when_four_occurs_first():
    assert canonicalize_residual_coloring([0, 4, 1, 3, 4]) == (0, 3, 1, 4, 3)
    assert canonicalize_residual_coloring([0, 3, 1, 4]) == (0, 3, 1, 4)


def test_true_variables_include_exact_colors_and_seen_color3_prefix():
    # Three vertices, five colors, then three residual-prefix auxiliaries.
    true = phase_true_variables([0, 3, 4], colors=5)

    assert true == {1, 9, 15, 17, 18}


def test_rewrite_flips_literal_signs_and_preserves_header(tmp_path: Path):
    source = tmp_path / "source.cnf"
    target = tmp_path / "target.cnf"
    source.write_text("p cnf 3 2\n1 -2 0\n-1 3 0\n", encoding="ascii")

    stats = rewrite_cnf_for_false_phase(source, target, {1, 3})

    assert target.read_text(encoding="ascii") == (
        "p cnf 3 2\n-1 -2 0\n1 -3 0\n"
    )
    assert stats == {"variables": 3, "clauses": 2, "flipped_variables": 2}


def test_unflip_model_is_an_involution_on_variable_signs():
    transformed = (-1, 2, -3, 4)

    assert unflip_model(transformed, {1, 3}) == (1, 2, 3, 4)
