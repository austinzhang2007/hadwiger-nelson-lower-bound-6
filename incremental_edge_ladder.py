"""Solve a hard coloring instance by adding target edges incrementally."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import random
from threading import Timer
import time
from typing import Callable, Iterable, Sequence

from pysat.solvers import Solver

from coloring_sat import ColoringSAT, Edge, parse_dimacs_edge
from d4_blocker_cegis import _write_json


def solve_with_timeout(
    solver: Solver,
    *,
    assumptions: Sequence[int],
    timeout_seconds: float,
) -> bool | None:
    """Run an interruptible limited solve, returning ``None`` on timeout."""

    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    timer = Timer(timeout_seconds, solver.interrupt)
    timer.start()
    try:
        return solver.solve_limited(
            assumptions=list(assumptions),
            expect_interrupt=True,
        )
    finally:
        timer.cancel()
        solver.clear_interrupt()


@dataclass(frozen=True, slots=True)
class EdgeLadderResult:
    status: str
    coloring: tuple[int, ...] | None
    active_edges: int
    target_edges: int
    sat_calls: int
    elapsed_seconds: float
    iterations: tuple[dict[str, object], ...]


def _edge_clauses(edge: Edge, color_count: int) -> list[list[int]]:
    left, right = edge
    return [
        [-(left * color_count + color + 1), -(right * color_count + color + 1)]
        for color in range(color_count)
    ]


def _decode_model(
    model: Sequence[int], vertex_count: int, color_count: int
) -> tuple[int, ...]:
    positive = {literal for literal in model if literal > 0}
    coloring = []
    for vertex in range(vertex_count):
        selected = [
            color
            for color in range(color_count)
            if vertex * color_count + color + 1 in positive
        ]
        if len(selected) != 1:
            raise RuntimeError("SAT model does not select exactly one color")
        coloring.append(selected[0])
    return tuple(coloring)


def solve_incremental_edge_ladder(
    *,
    vertex_count: int,
    base_edges: Iterable[Edge],
    target_edges: Iterable[Edge],
    color_count: int,
    initial_coloring: Sequence[int],
    seed: int,
    symmetry_clique: Sequence[int] = (),
    solver_name: str = "cadical195",
    max_satisfied_per_round: int | None = None,
    resume_active_edges: Iterable[Edge] | None = None,
    initial_sat_calls: int = 0,
    checkpoint_every: int | None = None,
    checkpoint: Callable[[set[Edge], tuple[int, ...], int], None] | None = None,
    solve_timeout_seconds: float | None = None,
    max_timed_out_edges: int = 8,
    progress: Callable[[str, dict[str, object]], None] | None = None,
) -> EdgeLadderResult:
    """Reach the target graph while preserving a nearby SAT model when possible.

    Every round adds a bounded number of not-yet-active target edges already
    satisfied by the current coloring and one randomly selected conflicting
    edge.  A bound of zero enables conflict-only rounds; once one model satisfies
    every remaining target edge, those edges are activated without another solve.
    """

    started = time.monotonic()
    if max_satisfied_per_round is not None and max_satisfied_per_round < 0:
        raise ValueError("max_satisfied_per_round must be nonnegative")
    if checkpoint_every is not None and checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive")
    if initial_sat_calls < 0:
        raise ValueError("initial_sat_calls must be nonnegative")
    if solve_timeout_seconds is not None:
        if solver_name != "glucose4":
            raise ValueError("timed edge trials currently require glucose4")
        if max_satisfied_per_round != 0:
            raise ValueError("timed edge trials require conflict-only rounds")
        if max_timed_out_edges < 1:
            raise ValueError("max_timed_out_edges must be positive")
    normalized_base = {
        (min(int(left), int(right)), max(int(left), int(right)))
        for left, right in base_edges
    }
    normalized_target = {
        (min(int(left), int(right)), max(int(left), int(right)))
        for left, right in target_edges
    }
    if not normalized_base <= normalized_target:
        raise ValueError("base edges must be a subset of target edges")
    if resume_active_edges is None:
        active = set(normalized_base)
    else:
        active = {
            (min(int(left), int(right)), max(int(left), int(right)))
            for left, right in resume_active_edges
        }
        if not normalized_base <= active <= normalized_target:
            raise ValueError(
                "resume active edges must lie between base and target edges"
            )
    encoding = ColoringSAT(
        vertex_count,
        active,
        color_count,
        break_color_symmetry=False,
        symmetry_clique=symmetry_clique,
        solver_name=solver_name,
    )
    coloring = tuple(int(color) for color in initial_coloring)
    if not encoding.validate(coloring):
        raise ValueError("initial coloring is invalid on the active graph")

    rng = random.Random(seed)
    remaining = normalized_target - active
    history: list[dict[str, object]] = []
    sat_calls = initial_sat_calls
    last_checkpoint_sat_calls: int | None = None
    timed_out_edges: set[Edge] = set()
    next_selector = encoding.variable_count + 1
    solver = Solver(name=solver_name, bootstrap_with=encoding.clauses)
    try:
        while remaining:
            all_satisfied = {
                edge for edge in remaining if coloring[edge[0]] != coloring[edge[1]]
            }
            if (
                max_satisfied_per_round is not None
                and len(all_satisfied) > max_satisfied_per_round
            ):
                satisfied = set(
                    rng.sample(
                        sorted(all_satisfied),
                        max_satisfied_per_round,
                    )
                )
            else:
                satisfied = all_satisfied
            if satisfied:
                solver.append_formula(
                    [
                        clause
                        for edge in sorted(satisfied)
                        for clause in _edge_clauses(edge, color_count)
                    ]
                )
                active.update(satisfied)
                remaining.difference_update(satisfied)
            if not remaining:
                break

            all_conflicting = sorted(
                edge
                for edge in remaining
                if coloring[edge[0]] == coloring[edge[1]]
            )
            if not all_conflicting:
                active.update(remaining)
                remaining.clear()
                break
            conflicting = [
                edge for edge in all_conflicting if edge not in timed_out_edges
            ]
            if not conflicting:
                return EdgeLadderResult(
                    status="UNKNOWN_TIMEOUT",
                    coloring=coloring,
                    active_edges=len(active),
                    target_edges=len(normalized_target),
                    sat_calls=sat_calls,
                    elapsed_seconds=time.monotonic() - started,
                    iterations=tuple(history),
                )
            chosen = conflicting[rng.randrange(len(conflicting))]
            phases = [
                vertex * color_count + color + 1
                if color == coloring[vertex]
                else -(vertex * color_count + color + 1)
                for vertex in range(vertex_count)
                for color in range(color_count)
            ]
            solver.set_phases(phases)
            solve_started = time.monotonic()
            selector: int | None = None
            if solve_timeout_seconds is None:
                solver.append_formula(_edge_clauses(chosen, color_count))
                satisfiable = solver.solve()
            else:
                selector = next_selector
                next_selector += 1
                solver.append_formula(
                    [
                        [-selector, *clause]
                        for clause in _edge_clauses(chosen, color_count)
                    ]
                )
                satisfiable = solve_with_timeout(
                    solver,
                    assumptions=(selector,),
                    timeout_seconds=solve_timeout_seconds,
                )
            solve_seconds = time.monotonic() - solve_started
            sat_calls += 1
            if satisfiable is None:
                if selector is None:
                    raise AssertionError("untimed solve returned unknown")
                solver.add_clause([-selector])
                timed_out_edges.add(chosen)
                record = {
                    "sat_call": sat_calls,
                    "satisfied_edges_added": len(satisfied),
                    "satisfied_edges_available": len(all_satisfied),
                    "conflicting_edges_before_choice": len(all_conflicting),
                    "chosen_conflicting_edge": list(chosen),
                    "active_edges": len(active),
                    "remaining_edges": len(remaining),
                    "solve_seconds": solve_seconds,
                    "satisfiable": None,
                    "timed_out": True,
                }
                history.append(record)
                if progress is not None:
                    progress("sat_call", record)
                if len(timed_out_edges) >= max_timed_out_edges:
                    return EdgeLadderResult(
                        status="UNKNOWN_TIMEOUT",
                        coloring=coloring,
                        active_edges=len(active),
                        target_edges=len(normalized_target),
                        sat_calls=sat_calls,
                        elapsed_seconds=time.monotonic() - started,
                        iterations=tuple(history),
                    )
                solver.delete()
                encoding = ColoringSAT(
                    vertex_count,
                    active,
                    color_count,
                    break_color_symmetry=False,
                    symmetry_clique=symmetry_clique,
                    solver_name=solver_name,
                )
                solver = Solver(name=solver_name, bootstrap_with=encoding.clauses)
                next_selector = encoding.variable_count + 1
                continue
            if selector is not None:
                solver.add_clause([selector])
            active.add(chosen)
            remaining.remove(chosen)
            record: dict[str, object] = {
                "sat_call": sat_calls,
                "satisfied_edges_added": len(satisfied),
                "satisfied_edges_available": len(all_satisfied),
                "conflicting_edges_before_choice": len(all_conflicting),
                "chosen_conflicting_edge": list(chosen),
                "active_edges": len(active),
                "remaining_edges": len(remaining),
                "solve_seconds": solve_seconds,
                "satisfiable": satisfiable,
                "timed_out": False,
            }
            history.append(record)
            if progress is not None:
                progress("sat_call", record)
            if not satisfiable:
                return EdgeLadderResult(
                    status="UNSAT_WITHOUT_PROOF",
                    coloring=None,
                    active_edges=len(active),
                    target_edges=len(normalized_target),
                    sat_calls=sat_calls,
                    elapsed_seconds=time.monotonic() - started,
                    iterations=tuple(history),
                )
            model = solver.get_model()
            if model is None:
                raise RuntimeError("SAT solver returned no model")
            coloring = _decode_model(model, vertex_count, color_count)
            if any(coloring[left] == coloring[right] for left, right in active):
                raise RuntimeError("solver model violates an active edge")
            timed_out_edges.clear()
            if (
                checkpoint is not None
                and (
                    checkpoint_every is None
                    or sat_calls % checkpoint_every == 0
                )
            ):
                checkpoint(set(active), coloring, sat_calls)
                last_checkpoint_sat_calls = sat_calls
    finally:
        solver.delete()

    if any(coloring[left] == coloring[right] for left, right in normalized_target):
        raise RuntimeError("final coloring violates a target edge")
    if checkpoint is not None and last_checkpoint_sat_calls != sat_calls:
        checkpoint(set(active), coloring, sat_calls)
    return EdgeLadderResult(
        status="SAT",
        coloring=coloring,
        active_edges=len(active),
        target_edges=len(normalized_target),
        sat_calls=sat_calls,
        elapsed_seconds=time.monotonic() - started,
        iterations=tuple(history),
    )


def solve_report_edge_ladder(
    *,
    base_report_path: str | Path,
    target_report_path: str | Path,
    output_path: str | Path,
    seed: int,
    solver_name: str = "cadical195",
    max_satisfied_per_round: int | None = None,
    checkpoint_path: str | Path | None = None,
    checkpoint_every: int = 100,
    resume_checkpoint_path: str | Path | None = None,
    solve_timeout_seconds: float | None = None,
    max_timed_out_edges: int = 8,
    progress: Callable[[str, dict[str, object]], None] | None = None,
) -> dict[str, object]:
    """Run the ladder between two same-vertex checkpoint reports."""

    with Path(base_report_path).open(encoding="utf-8") as handle:
        base_report = json.load(handle)
    coloring = base_report.get("whole_graph_coloring")
    if not isinstance(coloring, list):
        raise ValueError("base report has no coloring")
    base_edge_path = str(base_report["edge_path"])
    del base_report

    with Path(target_report_path).open(encoding="utf-8") as handle:
        target_report = json.load(handle)
    target_edge_path = str(target_report["edge_path"])
    symmetry_clique = tuple(
        int(vertex) for vertex in target_report.get("symmetry_clique_zero_based", ())
    )
    del target_report

    base_vertices, base_edges = parse_dimacs_edge(base_edge_path)
    target_vertices, target_edges = parse_dimacs_edge(target_edge_path)
    if base_vertices != target_vertices:
        raise ValueError("edge ladder reports must have the same vertex set")
    resume_active_edges: set[Edge] | None = None
    initial_sat_calls = 0
    if resume_checkpoint_path is not None:
        with Path(resume_checkpoint_path).open(encoding="utf-8") as handle:
            saved = json.load(handle)
        if saved.get("schema") != "edge-ladder-checkpoint-v1":
            raise ValueError("unsupported edge ladder checkpoint schema")
        if int(saved["vertex_count"]) != target_vertices:
            raise ValueError("checkpoint vertex count disagrees with target graph")
        if int(saved["color_count"]) != 5:
            raise ValueError("checkpoint color count is not five")
        raw_active = saved.get("active_edges")
        raw_coloring = saved.get("coloring")
        if not isinstance(raw_active, list) or not isinstance(raw_coloring, list):
            raise ValueError("checkpoint lacks active edges or coloring")
        resume_active_edges = {
            (min(int(edge[0]), int(edge[1])), max(int(edge[0]), int(edge[1])))
            for edge in raw_active
        }
        if len(resume_active_edges) != int(saved["active_edge_count"]):
            raise ValueError("checkpoint active edge count is inconsistent")
        coloring = [int(color) for color in raw_coloring]
        initial_sat_calls = int(saved["sat_calls"])

    checkpoint_callback = None
    if checkpoint_path is not None:
        checkpoint_path = Path(checkpoint_path)

        def persist_checkpoint(
            active_edges: set[Edge],
            active_coloring: tuple[int, ...],
            sat_calls: int,
        ) -> None:
            checkpoint_payload: dict[str, object] = {
                "schema": "edge-ladder-checkpoint-v1",
                "base_report": str(base_report_path),
                "target_report": str(target_report_path),
                "vertex_count": target_vertices,
                "color_count": 5,
                "active_edge_count": len(active_edges),
                "active_edges": [list(edge) for edge in sorted(active_edges)],
                "coloring": list(active_coloring),
                "sat_calls": sat_calls,
                "seed": seed,
                "solver": solver_name,
                "max_satisfied_per_round": max_satisfied_per_round,
                "solve_timeout_seconds": solve_timeout_seconds,
                "max_timed_out_edges": max_timed_out_edges,
            }
            _write_json(checkpoint_path, checkpoint_payload)
            if progress is not None:
                progress(
                    "checkpoint_written",
                    {
                        "checkpoint_path": str(checkpoint_path),
                        "active_edges": len(active_edges),
                        "sat_calls": sat_calls,
                    },
                )

        checkpoint_callback = persist_checkpoint
    result = solve_incremental_edge_ladder(
        vertex_count=base_vertices,
        base_edges=base_edges,
        target_edges=target_edges,
        color_count=5,
        initial_coloring=coloring,
        seed=seed,
        symmetry_clique=symmetry_clique,
        solver_name=solver_name,
        max_satisfied_per_round=max_satisfied_per_round,
        resume_active_edges=resume_active_edges,
        initial_sat_calls=initial_sat_calls,
        checkpoint_every=checkpoint_every,
        checkpoint=checkpoint_callback,
        solve_timeout_seconds=solve_timeout_seconds,
        max_timed_out_edges=max_timed_out_edges,
        progress=progress,
    )
    payload: dict[str, object] = {
        **asdict(result),
        "base_report": str(base_report_path),
        "target_report": str(target_report_path),
        "seed": seed,
        "solver": solver_name,
        "max_satisfied_per_round": max_satisfied_per_round,
        "checkpoint_path": str(checkpoint_path) if checkpoint_path else None,
        "checkpoint_every": checkpoint_every,
        "solve_timeout_seconds": solve_timeout_seconds,
        "max_timed_out_edges": max_timed_out_edges,
        "resumed_from_checkpoint": (
            str(resume_checkpoint_path) if resume_checkpoint_path else None
        ),
        "color_class_sizes": (
            [result.coloring.count(color) for color in range(5)]
            if result.coloring is not None
            else None
        ),
    }
    _write_json(Path(output_path), payload)
    return payload


def apply_ladder_result_to_report(
    *,
    ladder_result_path: str | Path,
    target_report_path: str | Path,
) -> dict[str, object]:
    """Validate a completed SAT ladder and attach its model to the target report."""

    ladder_result_path = Path(ladder_result_path)
    target_report_path = Path(target_report_path)
    with ladder_result_path.open(encoding="utf-8") as handle:
        ladder = json.load(handle)
    if ladder.get("status") != "SAT":
        raise ValueError("only a completed SAT ladder can update a target report")
    coloring = ladder.get("coloring")
    if not isinstance(coloring, list):
        raise ValueError("SAT ladder result has no coloring")
    if int(ladder["active_edges"]) != int(ladder["target_edges"]):
        raise ValueError("SAT ladder did not activate every target edge")

    with target_report_path.open(encoding="utf-8") as handle:
        target = json.load(handle)
    vertex_count, edges = parse_dimacs_edge(str(target["edge_path"]))
    if vertex_count != int(target["vertices"]):
        raise ValueError("target report vertex count disagrees with edge file")
    if len(edges) != int(target["edges"]):
        raise ValueError("target report edge count disagrees with edge file")
    if len(edges) != int(ladder["target_edges"]):
        raise ValueError("ladder target edge count disagrees with target report")
    normalized_coloring = tuple(int(color) for color in coloring)
    if not ColoringSAT(
        vertex_count,
        edges,
        5,
        break_color_symmetry=False,
    ).validate(normalized_coloring):
        raise ValueError("ladder coloring violates the target edge file")

    target["status"] = "EXACT_INTERACTIONS_AUGMENTED_WITH_5_COLORING"
    target["lower_bound_6_proved"] = False
    target["whole_graph_5_sat"] = True
    target["whole_graph_coloring"] = list(normalized_coloring)
    metadata_keys = (
        "status",
        "active_edges",
        "target_edges",
        "sat_calls",
        "elapsed_seconds",
        "seed",
        "solver",
        "max_satisfied_per_round",
        "color_class_sizes",
        "base_report",
        "target_report",
    )
    target["incremental_edge_ladder"] = {
        "result_path": str(ladder_result_path),
        **{key: ladder[key] for key in metadata_keys if key in ladder},
    }
    _write_json(target_report_path, target)
    return target


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-report", required=True)
    parser.add_argument("--target-report", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--max-satisfied-per-round", type=int)
    parser.add_argument("--checkpoint")
    parser.add_argument("--checkpoint-every", type=int, default=100)
    parser.add_argument("--resume-checkpoint")
    parser.add_argument("--solve-timeout-seconds", type=float)
    parser.add_argument("--max-timed-out-edges", type=int, default=8)
    args = parser.parse_args()

    def report(event: str, details: dict[str, object]) -> None:
        print(json.dumps({"event": event, **details}, sort_keys=True), flush=True)

    payload = solve_report_edge_ladder(
        base_report_path=args.base_report,
        target_report_path=args.target_report,
        output_path=args.output,
        seed=args.seed,
        solver_name=args.solver,
        max_satisfied_per_round=args.max_satisfied_per_round,
        checkpoint_path=args.checkpoint,
        checkpoint_every=args.checkpoint_every,
        resume_checkpoint_path=args.resume_checkpoint,
        solve_timeout_seconds=args.solve_timeout_seconds,
        max_timed_out_edges=args.max_timed_out_edges,
        progress=report,
    )
    print(
        json.dumps(
            {
                "status": payload["status"],
                "active_edges": payload["active_edges"],
                "target_edges": payload["target_edges"],
                "sat_calls": payload["sat_calls"],
                "elapsed_seconds": payload["elapsed_seconds"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    _main()
