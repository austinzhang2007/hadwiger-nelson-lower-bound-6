import json
from pathlib import Path

from graph_verifier import verify_drat
import pytest
from pysat.solvers import pysolvers

from proof_solver import (
    decode_graph_coloring_model,
    export_validated_coloring,
    solve_cnf,
)


ROOT = Path(__file__).resolve().parents[1]


def test_proof_solver_saves_a_sat_model(tmp_path: Path) -> None:
    cnf = tmp_path / "sat.cnf"
    cnf.write_text("p cnf 1 1\n1 0\n", encoding="ascii")
    result_path = tmp_path / "sat.json"

    result = solve_cnf(
        cnf_path=cnf,
        result_path=result_path,
        proof_path=None,
        solver_name="cadical195",
    )

    assert result["status"] == "SAT"
    assert result["model"] == [1]
    with result_path.open(encoding="utf-8") as handle:
        assert json.load(handle)["status"] == "SAT"


def test_proof_solver_saves_a_drat_trim_verified_unsat_proof(
    tmp_path: Path,
) -> None:
    cnf = tmp_path / "unsat.cnf"
    cnf.write_text("p cnf 1 2\n1 0\n-1 0\n", encoding="ascii")
    result_path = tmp_path / "unsat.json"
    proof_path = tmp_path / "unsat.drat"

    result = solve_cnf(
        cnf_path=cnf,
        result_path=result_path,
        proof_path=proof_path,
        solver_name="glucose4",
    )
    proof = verify_drat(
        ROOT / "third_party/drat-trim/drat-trim",
        cnf,
        proof_path,
    )

    assert result["status"] == "UNSAT_UNVERIFIED"
    assert proof_path.is_file()
    assert proof.valid


def test_decode_graph_coloring_model_validates_every_edge() -> None:
    model = [1, -2, -3, -4, 5, -6]

    coloring = decode_graph_coloring_model(
        model,
        vertex_count=2,
        color_count=3,
        edges={(0, 1)},
    )

    assert coloring == (0, 1)


def test_decode_graph_coloring_model_rejects_a_conflict() -> None:
    with pytest.raises(ValueError, match="monochromatic"):
        decode_graph_coloring_model(
            [1, -2, -3, 4, -5, -6],
            vertex_count=2,
            color_count=3,
            edges={(0, 1)},
        )


def test_proof_solver_records_pysat_wrapped_interrupt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class InterruptingSolver:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def solve(self):
            raise pysolvers.error("Caught keyboard interrupt")

    cnf = tmp_path / "unknown.cnf"
    cnf.write_text("p cnf 1 1\n1 0\n", encoding="ascii")
    result_path = tmp_path / "unknown.json"
    monkeypatch.setattr(
        "proof_solver.Solver",
        lambda **kwargs: InterruptingSolver(),
    )

    result = solve_cnf(
        cnf_path=cnf,
        result_path=result_path,
        proof_path=None,
        solver_name="cadical195",
    )

    assert result["status"] == "UNKNOWN_INTERRUPTED"


def test_export_validated_coloring_writes_coloring_and_phases(
    tmp_path: Path,
) -> None:
    result_path = tmp_path / "result.json"
    result_path.write_text(
        json.dumps(
            {
                "status": "SAT",
                "solver": "cadical195",
                "model": [1, -2, -3, -4, 5, -6],
            }
        ),
        encoding="utf-8",
    )
    coloring_path = tmp_path / "coloring.json"
    phase_path = tmp_path / "phases.json"

    payload = export_validated_coloring(
        result_path=result_path,
        coloring_path=coloring_path,
        phase_path=phase_path,
        vertex_count=2,
        color_count=3,
        edges={(0, 1)},
    )

    assert payload["coloring"] == [0, 1]
    assert payload["class_sizes"] == [1, 1, 0]
    with phase_path.open(encoding="utf-8") as handle:
        assert json.load(handle)["phases"] == [1, -2, -3, -4, 5, -6]
