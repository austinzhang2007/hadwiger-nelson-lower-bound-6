import json
import time

from pysat.examples.genhard import PHP
from pysat.solvers import Solver

from coloring_sat import ColoringSAT, write_dimacs_edge
from incremental_edge_ladder import (
    apply_ladder_result_to_report,
    solve_with_timeout,
    solve_incremental_edge_ladder,
    solve_report_edge_ladder,
)


def test_incremental_edge_ladder_reaches_full_sat_graph() -> None:
    target = {(0, 1), (1, 2), (0, 2)}

    result = solve_incremental_edge_ladder(
        vertex_count=3,
        base_edges=set(),
        target_edges=target,
        color_count=3,
        initial_coloring=(0, 0, 0),
        seed=17,
    )

    assert result.status == "SAT"
    assert result.active_edges == 3
    assert result.coloring is not None
    assert ColoringSAT(3, target, 3, break_color_symmetry=False).validate(
        result.coloring
    )


def test_incremental_edge_ladder_detects_unsat_prefix() -> None:
    target = {
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 3),
    }

    result = solve_incremental_edge_ladder(
        vertex_count=4,
        base_edges=set(),
        target_edges=target,
        color_count=3,
        initial_coloring=(0, 0, 0, 0),
        seed=23,
    )

    assert result.status == "UNSAT_WITHOUT_PROOF"
    assert result.coloring is None
    assert result.active_edges <= len(target)


def test_apply_sat_ladder_result_to_target_report(tmp_path) -> None:
    edge_path = tmp_path / "target.edge"
    target_path = tmp_path / "target.json"
    ladder_path = tmp_path / "ladder.json"
    edges = {(0, 1), (1, 2), (0, 2)}
    write_dimacs_edge(edge_path, 3, edges)
    target_path.write_text(
        json.dumps(
            {
                "status": "EXACT_INTERACTIONS_PENDING_SAT",
                "vertices": 3,
                "edges": 3,
                "edge_path": str(edge_path),
                "whole_graph_5_sat": False,
                "whole_graph_coloring": None,
            }
        ),
        encoding="utf-8",
    )
    ladder_path.write_text(
        json.dumps(
            {
                "status": "SAT",
                "coloring": [0, 1, 2],
                "active_edges": 3,
                "target_edges": 3,
                "sat_calls": 2,
                "elapsed_seconds": 0.25,
                "seed": 17,
                "solver": "cadical195",
                "max_satisfied_per_round": 1,
            }
        ),
        encoding="utf-8",
    )

    updated = apply_ladder_result_to_report(
        ladder_result_path=ladder_path,
        target_report_path=target_path,
    )

    persisted = json.loads(target_path.read_text(encoding="utf-8"))
    assert updated == persisted
    assert persisted["status"] == "EXACT_INTERACTIONS_AUGMENTED_WITH_5_COLORING"
    assert persisted["whole_graph_5_sat"] is True
    assert persisted["whole_graph_coloring"] == [0, 1, 2]
    assert persisted["incremental_edge_ladder"]["sat_calls"] == 2


def test_incremental_edge_ladder_resumes_from_sat_checkpoint() -> None:
    target = {(0, 1), (1, 2), (0, 2)}
    snapshots: list[tuple[set[tuple[int, int]], tuple[int, ...], int]] = []

    solve_incremental_edge_ladder(
        vertex_count=3,
        base_edges=set(),
        target_edges=target,
        color_count=3,
        initial_coloring=(0, 0, 0),
        seed=31,
        max_satisfied_per_round=1,
        checkpoint_every=1,
        checkpoint=lambda edges, coloring, calls: snapshots.append(
            (set(edges), tuple(coloring), calls)
        ),
    )
    active, coloring, sat_calls = snapshots[0]

    resumed = solve_incremental_edge_ladder(
        vertex_count=3,
        base_edges=set(),
        target_edges=target,
        color_count=3,
        initial_coloring=coloring,
        seed=32,
        max_satisfied_per_round=1,
        resume_active_edges=active,
        initial_sat_calls=sat_calls,
    )

    assert resumed.status == "SAT"
    assert resumed.active_edges == len(target)
    assert resumed.sat_calls >= sat_calls
    assert resumed.coloring is not None
    assert ColoringSAT(3, target, 3, break_color_symmetry=False).validate(
        resumed.coloring
    )


def test_report_ladder_persists_and_resumes_checkpoint(tmp_path) -> None:
    base_edge_path = tmp_path / "base.edge"
    target_edge_path = tmp_path / "target.edge"
    base_report_path = tmp_path / "base.json"
    target_report_path = tmp_path / "target.json"
    checkpoint_path = tmp_path / "checkpoint.json"
    first_output_path = tmp_path / "first.json"
    resumed_output_path = tmp_path / "resumed.json"
    target = {(0, 1), (1, 2), (0, 2)}
    write_dimacs_edge(base_edge_path, 3, set())
    write_dimacs_edge(target_edge_path, 3, target)
    base_report_path.write_text(
        json.dumps(
            {
                "edge_path": str(base_edge_path),
                "whole_graph_coloring": [0, 0, 0],
            }
        ),
        encoding="utf-8",
    )
    target_report_path.write_text(
        json.dumps(
            {
                "edge_path": str(target_edge_path),
                "symmetry_clique_zero_based": [],
            }
        ),
        encoding="utf-8",
    )

    first = solve_report_edge_ladder(
        base_report_path=base_report_path,
        target_report_path=target_report_path,
        output_path=first_output_path,
        seed=41,
        max_satisfied_per_round=1,
        checkpoint_path=checkpoint_path,
        checkpoint_every=1,
    )
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert first["status"] == "SAT"
    assert checkpoint["schema"] == "edge-ladder-checkpoint-v1"
    assert checkpoint["active_edge_count"] == len(target)

    resumed = solve_report_edge_ladder(
        base_report_path=base_report_path,
        target_report_path=target_report_path,
        output_path=resumed_output_path,
        seed=42,
        max_satisfied_per_round=1,
        checkpoint_path=checkpoint_path,
        checkpoint_every=1,
        resume_checkpoint_path=checkpoint_path,
    )

    assert resumed["status"] == "SAT"
    assert resumed["active_edges"] == len(target)
    assert resumed["sat_calls"] >= checkpoint["sat_calls"]

    timed = solve_report_edge_ladder(
        base_report_path=base_report_path,
        target_report_path=target_report_path,
        output_path=tmp_path / "timed.json",
        seed=43,
        solver_name="glucose4",
        max_satisfied_per_round=0,
        solve_timeout_seconds=1.0,
        max_timed_out_edges=2,
    )
    assert timed["status"] == "SAT"


def test_incremental_edge_ladder_supports_conflict_only_rounds() -> None:
    target = {(0, 1), (1, 2), (0, 2)}

    result = solve_incremental_edge_ladder(
        vertex_count=3,
        base_edges=set(),
        target_edges=target,
        color_count=3,
        initial_coloring=(0, 0, 0),
        seed=53,
        max_satisfied_per_round=0,
    )

    assert result.status == "SAT"
    assert result.active_edges == len(target)
    assert result.coloring is not None
    assert all(
        iteration["satisfied_edges_added"] == 0
        for iteration in result.iterations
    )
    assert ColoringSAT(3, target, 3, break_color_symmetry=False).validate(
        result.coloring
    )


def test_glucose_solve_timeout_returns_unknown() -> None:
    started = time.monotonic()
    with Solver(name="glucose4", bootstrap_with=PHP(nof_holes=20).clauses) as solver:
        result = solve_with_timeout(
            solver,
            assumptions=(),
            timeout_seconds=0.01,
        )

    assert result is None
    assert time.monotonic() - started < 2.0


def test_ladder_returns_unknown_without_activating_timed_out_edge(
    monkeypatch,
) -> None:
    created = []

    class FakeSolver:
        def __init__(self, *args, **kwargs) -> None:
            self.deleted = False
            created.append(self)

        def __enter__(self):
            return self

        def __exit__(self, *args) -> None:
            self.delete()

        def append_formula(self, clauses) -> None:
            pass

        def set_phases(self, phases) -> None:
            pass

        def add_clause(self, clause) -> None:
            pass

        def delete(self) -> None:
            self.deleted = True

    monkeypatch.setattr("incremental_edge_ladder.Solver", FakeSolver)
    monkeypatch.setattr(
        "incremental_edge_ladder.solve_with_timeout",
        lambda solver, *, assumptions, timeout_seconds: None,
    )

    result = solve_incremental_edge_ladder(
        vertex_count=3,
        base_edges=set(),
        target_edges={(0, 1), (1, 2), (0, 2)},
        color_count=3,
        initial_coloring=(0, 0, 0),
        seed=59,
        solver_name="glucose4",
        max_satisfied_per_round=0,
        solve_timeout_seconds=0.01,
        max_timed_out_edges=2,
    )

    assert result.status == "UNKNOWN_TIMEOUT"
    assert result.active_edges == 0
    assert result.coloring == (0, 0, 0)
    assert result.iterations[-1]["timed_out"] is True
    assert len(created) == 2
    assert all(solver.deleted for solver in created)
