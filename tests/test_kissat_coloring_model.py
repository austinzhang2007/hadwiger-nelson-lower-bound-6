from pathlib import Path

import pytest

from kissat_coloring_model import (
    parse_kissat_output,
    validate_kissat_lean_coloring,
)


def test_parse_sat_model_across_multiple_value_lines(tmp_path: Path):
    log = tmp_path / "kissat.log"
    log.write_text(
        "c example\ns SATISFIABLE\nv 1 2 -3\nv -4 -5 6 0\n",
        encoding="ascii",
    )

    status, model = parse_kissat_output(log)

    assert status == "SAT"
    assert model == (1, 2, -3, -4, -5, 6)


def test_parse_unsat_has_no_model(tmp_path: Path):
    log = tmp_path / "kissat.log"
    log.write_text("s UNSATISFIABLE\n", encoding="ascii")

    assert parse_kissat_output(log) == ("UNSAT", ())


def test_validate_sat_lean_model_rechecks_all_edges(tmp_path: Path):
    log = tmp_path / "kissat.log"
    log.write_text(
        "s SATISFIABLE\nv 1 2 -3 -4 -5 6 0\n",
        encoding="ascii",
    )

    coloring = validate_kissat_lean_coloring(
        log,
        vertex_count=2,
        edges={(0, 1)},
        colors=3,
    )

    assert coloring == (0, 2)


def test_incomplete_running_log_is_rejected(tmp_path: Path):
    log = tmp_path / "kissat.log"
    log.write_text("c still running\n", encoding="ascii")

    with pytest.raises(ValueError, match="terminal status"):
        parse_kissat_output(log)
