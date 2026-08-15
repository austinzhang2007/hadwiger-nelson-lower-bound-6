"""Second-layer unit-circle closure CEGIS."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp

from cegis import CertifiedCandidate, generate_certified_unit_circle_candidates
from coloring_sat import ColoringSAT
from configuration_cegis import (
    _candidate_from_report,
    augmented_edges_from_pool,
    certified_candidate_unit_edges,
    configuration_cegis_step,
)
from exact_geometry import Point, complete_unit_edges, serialize_graph
from graph_verifier import load_graph_json


def write_checkpoint(path: str | Path, report: dict[str, object]) -> None:
    """Atomically replace a JSON checkpoint."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(destination)


def load_exact_colored_graph(
    graph_path: str | Path,
    report_path: str | Path,
    *,
    color_count: int,
) -> tuple[list[Point], set[tuple[int, int]], tuple[int, ...]]:
    """Load an exact graph and validate the saved whole-pool coloring."""

    points, edges = load_graph_json(graph_path)
    with Path(report_path).open(encoding="utf-8") as handle:
        report = json.load(handle)
    coloring = tuple(
        int(color)
        for color in report["whole_pool_failure_certificate"]["coloring"]
    )
    if not ColoringSAT(
        len(points),
        edges,
        color_count,
        break_color_symmetry=False,
    ).validate(coloring):
        raise ValueError("saved coloring is invalid on the exact base graph")
    return points, edges, coloring


def solve_full_candidate_pool(
    *,
    base_vertex_count: int,
    base_edges: Iterable[tuple[int, int]],
    candidates: Sequence[CertifiedCandidate],
    candidate_edges: Iterable[tuple[int, int]],
    color_count: int,
) -> tuple[set[tuple[int, int]], tuple[int, ...] | None]:
    """Solve the entire augmented pool without fixing the base coloring."""

    _, edges = augmented_edges_from_pool(
        base_vertex_count,
        base_edges,
        candidates,
        candidate_edges,
        range(len(candidates)),
    )
    coloring = ColoringSAT(
        base_vertex_count + len(candidates),
        edges,
        color_count,
    ).solve()
    return edges, coloring


def finalize_closure_checkpoint(
    graph_path: str | Path,
    checkpoint_path: str | Path,
    *,
    color_count: int = 5,
) -> dict[str, object]:
    """Solve and serialize the full listed candidate pool from a checkpoint."""

    started = time.monotonic()
    graph_path = Path(graph_path)
    checkpoint_path = Path(checkpoint_path)
    points, base_edges = load_graph_json(graph_path)
    with checkpoint_path.open(encoding="utf-8") as handle:
        report = json.load(handle)
    candidates = [
        _candidate_from_report(record) for record in report["candidate_pool"]
    ]
    candidate_edges = {
        tuple(int(index) for index in edge)
        for edge in report["candidate_unit_edges_zero_based"]
    }
    listed_edges, coloring = solve_full_candidate_pool(
        base_vertex_count=len(points),
        base_edges=base_edges,
        candidates=candidates,
        candidate_edges=candidate_edges,
        color_count=color_count,
    )
    report["cegis_checkpoint_status_before_direct_solve"] = report["status"]
    if coloring is None:
        report["status"] = "FULL_LISTED_POOL_5_UNSAT_WITHOUT_PROOF"
        report["direct_full_pool_solve"] = {
            "sat": False,
            "elapsed_seconds": time.monotonic() - started,
        }
        write_checkpoint(checkpoint_path, report)
        return report

    full_points = [*points, *[candidate.point for candidate in candidates]]
    listed_graph_path = checkpoint_path.with_name(
        checkpoint_path.stem + "-listed.graph.json"
    )
    _write_graph(
        listed_graph_path,
        "Second unit-circle closure with certified listed edges",
        full_points,
        listed_edges,
    )
    report["status"] = "EXACT_LISTED_EDGE_WHOLE_POOL_5_COLORING"
    report["whole_pool_failure_certificate"] = {
        "listed_graph_path": str(listed_graph_path),
        "vertices": len(full_points),
        "certified_listed_edges": len(listed_edges),
        "coloring": list(coloring),
        "coloring_valid_on_listed_edges": True,
        "complete_edge_audit_attempted": False,
    }
    report["direct_full_pool_solve"] = {
        "sat": True,
        "elapsed_seconds": time.monotonic() - started,
    }
    write_checkpoint(checkpoint_path, report)
    return report


def run_direct_closure_audit(
    graph_path: str | Path,
    output_path: str | Path,
    *,
    minimum_candidate_degree: int = 3,
    color_count: int = 5,
    proof_checker_path: str | Path = "third_party/drat-trim/drat-trim",
) -> dict[str, object]:
    """Generate one closure, rebuild its complete graph, and solve it directly."""

    started = time.monotonic()
    graph_path = Path(graph_path)
    output_path = Path(output_path)
    points, _ = load_graph_json(graph_path)
    candidates, candidate_stats = generate_certified_unit_circle_candidates(
        points,
        min_neighbors=minimum_candidate_degree,
    )
    candidate_edges, candidate_edge_stats = certified_candidate_unit_edges(
        candidates
    )
    full_points = [*points, *[candidate.point for candidate in candidates]]
    if len(set(full_points)) != len(full_points):
        raise RuntimeError("direct closure contains duplicate exact points")
    complete_edges = complete_unit_edges(full_points)
    encoding = ColoringSAT(len(full_points), complete_edges, color_count)
    coloring = encoding.solve()
    graph_output_path = output_path.with_name(
        output_path.stem + "-complete.graph.json"
    )
    _write_graph(
        graph_output_path,
        "Direct exact unit-circle closure audit",
        full_points,
        complete_edges,
    )
    report: dict[str, object] = {
        "status": "RUNNING",
        "lower_bound_6_proved": False,
        "base_graph_path": str(graph_path),
        "minimum_candidate_degree": minimum_candidate_degree,
        "color_count": color_count,
        "candidate_generation": asdict(candidate_stats),
        "candidate_interactions": asdict(candidate_edge_stats),
        "candidate_pool": [
            _candidate_record(candidate, color_count=color_count)
            for candidate in candidates
        ],
        "complete_graph": {
            "path": str(graph_output_path),
            "vertices": len(full_points),
            "edges": len(complete_edges),
        },
    }
    if coloring is not None:
        report["status"] = "EXACT_COMPLETE_WHOLE_POOL_COLORING"
        report["whole_pool_failure_certificate"] = {
            "coloring": list(coloring),
            "coloring_valid_on_exact_complete_edges": True,
        }
    else:
        cnf_path = output_path.with_suffix(".cnf")
        proof_path = output_path.with_suffix(".drat")
        encoding.write_dimacs_cnf(cnf_path)
        proof_result = encoding.solve_with_drat(proof_path)
        if proof_result is not None:
            raise RuntimeError("proof solver disagrees with direct closure UNSAT")
        proof_verified, transcript = _verify_drat(
            Path(proof_checker_path),
            cnf_path,
            proof_path,
        )
        report["status"] = (
            "EXACT_COMPLETE_UNSAT_DRAT_VERIFIED"
            if proof_verified
            else "EXACT_COMPLETE_UNSAT_PROOF_VERIFICATION_FAILED"
        )
        report["potential_lower_bound_certificate"] = {
            "cnf_path": str(cnf_path),
            "proof_path": str(proof_path),
            "drat_verified": proof_verified,
            "drat_transcript": transcript,
        }
        report["lower_bound_6_proved"] = color_count == 5 and proof_verified
    report["elapsed_seconds"] = time.monotonic() - started
    write_checkpoint(output_path, report)
    return report


def _candidate_record(
    candidate: CertifiedCandidate,
    base_coloring: Sequence[int] | None = None,
    *,
    color_count: int = 5,
) -> dict[str, object]:
    x, y = candidate.point.approximate()
    record: dict[str, object] = {
        "x_exact": sp.sstr(candidate.point.x),
        "y_exact": sp.sstr(candidate.point.y),
        "x_approx": x,
        "y_approx": y,
        "neighbors_zero_based": list(candidate.neighbors),
        "degree_to_base": len(candidate.neighbors),
        "generator_pair_zero_based": list(candidate.generator_pair),
    }
    if base_coloring is not None:
        seen = sorted({int(base_coloring[index]) for index in candidate.neighbors})
        record["neighbor_colors"] = seen
        record["allowed_colors"] = sorted(set(range(color_count)).difference(seen))
    return record


def _write_graph(
    path: Path,
    name: str,
    points: Sequence[Point],
    edges: Iterable[tuple[int, int]],
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            serialize_graph(name, points, edges),
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")


def _verify_drat(
    checker: Path,
    cnf_path: Path,
    proof_path: Path,
) -> tuple[bool, str]:
    completed = subprocess.run(
        [str(checker), str(cnf_path), str(proof_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    transcript = completed.stdout + completed.stderr
    return completed.returncode == 0 and "s VERIFIED" in transcript, transcript


def run_closure_cegis(
    graph_path: str | Path,
    coloring_report_path: str | Path,
    output_path: str | Path,
    *,
    max_iterations: int = 100,
    minimum_candidate_degree: int = 3,
    color_count: int = 5,
    resume: bool = False,
    complete_audit: bool = False,
    proof_checker_path: str | Path = "third_party/drat-trim/drat-trim",
) -> dict[str, object]:
    """Run a checkpointed CEGIS layer over a saved exact colored graph."""

    started = time.monotonic()
    graph_path = Path(graph_path)
    coloring_report_path = Path(coloring_report_path)
    output_path = Path(output_path)
    points, base_edges, initial_coloring = load_exact_colored_graph(
        graph_path,
        coloring_report_path,
        color_count=color_count,
    )
    candidates, candidate_stats = generate_certified_unit_circle_candidates(
        points,
        min_neighbors=minimum_candidate_degree,
    )
    candidate_edges, candidate_edge_stats = certified_candidate_unit_edges(
        candidates
    )

    if resume and output_path.exists():
        with output_path.open(encoding="utf-8") as handle:
            report = json.load(handle)
        if report["base"]["vertices"] != len(points):
            raise ValueError("checkpoint base graph does not match requested graph")
        if report["candidate_generation"]["certified_candidates"] != len(candidates):
            raise ValueError("checkpoint candidate pool does not match regenerated pool")
        state = report["resume_state"]
        start_iteration = int(state["next_iteration"])
        selected = {int(index) for index in state["selected_candidate_indices"]}
        current_coloring = tuple(
            int(color) for color in state["current_base_coloring"]
        )
        previous_elapsed = float(report.get("elapsed_seconds", 0.0))
        report["status"] = "RUNNING"
    else:
        start_iteration = 0
        selected: set[int] = set()
        current_coloring = initial_coloring
        previous_elapsed = 0.0
        report: dict[str, object] = {
            "status": "RUNNING",
            "lower_bound_6_proved": False,
            "base": {
                "graph_path": str(graph_path),
                "coloring_report_path": str(coloring_report_path),
                "vertices": len(points),
                "exact_complete_edges": len(base_edges),
            },
            "minimum_candidate_degree": minimum_candidate_degree,
            "candidate_generation": asdict(candidate_stats),
            "candidate_interactions": asdict(candidate_edge_stats),
            "candidate_pool": [
                _candidate_record(candidate) for candidate in candidates
            ],
            "candidate_unit_edges_zero_based": [
                list(edge) for edge in sorted(candidate_edges)
            ],
            "initial_base_coloring": list(initial_coloring),
            "iterations": [],
        }

    termination = "ITERATION_LIMIT"
    for iteration_index in range(start_iteration, max_iterations):
        step_started = time.monotonic()
        fixed_coloring = current_coloring
        step = configuration_cegis_step(
            base_vertex_count=len(points),
            base_edges=base_edges,
            candidates=candidates,
            candidate_edges=candidate_edges,
            selected_indices=selected,
            base_coloring=fixed_coloring,
            color_count=color_count,
        )
        core_set = set(step.blocking_core)
        iteration: dict[str, object] = {
            "iteration": iteration_index,
            "result": step.result,
            "fixed_base_coloring": list(fixed_coloring),
            "blocking_core_candidate_indices": list(step.blocking_core),
            "blocking_core": [
                _candidate_record(
                    candidates[index],
                    fixed_coloring,
                    color_count=color_count,
                )
                for index in step.blocking_core
            ],
            "blocking_core_unit_edges": [
                list(edge)
                for edge in sorted(candidate_edges)
                if edge[0] in core_set and edge[1] in core_set
            ],
            "selected_candidate_indices": list(step.selected_indices),
            "selected_candidates": len(step.selected_indices),
            "augmented_edges": step.augmented_edges,
            "elapsed_seconds": time.monotonic() - step_started,
        }
        report["iterations"].append(iteration)

        if step.result == "SAT_COUNTEREXAMPLE":
            if step.new_base_coloring is None:
                raise RuntimeError("SAT counterexample omitted its base coloring")
            selected = set(step.selected_indices)
            current_coloring = step.new_base_coloring
            iteration["new_base_coloring"] = list(current_coloring)
            report["resume_state"] = {
                "next_iteration": iteration_index + 1,
                "selected_candidate_indices": sorted(selected),
                "current_base_coloring": list(current_coloring),
            }
            report["elapsed_seconds"] = (
                previous_elapsed + time.monotonic() - started
            )
            write_checkpoint(output_path, report)
            continue

        if step.result == "WHOLE_POOL_EXTENDS":
            if step.full_pool_assignment is None:
                raise RuntimeError("whole-pool extension omitted its assignment")
            full_order, listed_edges = augmented_edges_from_pool(
                len(points),
                base_edges,
                candidates,
                candidate_edges,
                range(len(candidates)),
            )
            full_points = [
                *points,
                *[candidates[index].point for index in full_order],
            ]
            full_coloring = (
                *fixed_coloring,
                *step.full_pool_assignment,
            )
            listed_coloring_valid = ColoringSAT(
                len(full_points),
                listed_edges,
                color_count,
                break_color_symmetry=False,
            ).validate(full_coloring)
            if not listed_coloring_valid:
                raise RuntimeError("whole-pool model violates a certified listed edge")
            listed_graph_path = output_path.with_name(
                output_path.stem + "-listed.graph.json"
            )
            _write_graph(
                listed_graph_path,
                "Second unit-circle closure with certified listed edges",
                full_points,
                listed_edges,
            )
            certificate: dict[str, object] = {
                "listed_graph_path": str(listed_graph_path),
                "vertices": len(full_points),
                "certified_listed_edges": len(listed_edges),
                "coloring": list(full_coloring),
                "coloring_valid_on_listed_edges": True,
                "complete_edge_audit_attempted": complete_audit,
            }
            termination = "EXACT_LISTED_EDGE_WHOLE_POOL_5_COLORING"
            if complete_audit:
                complete_edges = complete_unit_edges(full_points)
                complete_coloring_valid = ColoringSAT(
                    len(full_points),
                    complete_edges,
                    color_count,
                    break_color_symmetry=False,
                ).validate(full_coloring)
                complete_graph_path = output_path.with_name(
                    output_path.stem + "-complete.graph.json"
                )
                _write_graph(
                    complete_graph_path,
                    "Second unit-circle closure complete unit-distance graph",
                    full_points,
                    complete_edges,
                )
                certificate.update(
                    {
                        "complete_graph_path": str(complete_graph_path),
                        "exact_complete_edges": len(complete_edges),
                        "coloring_valid_on_exact_complete_edges": (
                            complete_coloring_valid
                        ),
                    }
                )
                termination = (
                    "EXACT_COMPLETE_WHOLE_POOL_5_COLORING"
                    if complete_coloring_valid
                    else "LISTED_COLORING_FAILED_COMPLETE_EDGE_AUDIT"
                )
            report["whole_pool_failure_certificate"] = certificate
            break

        if step.result != "AUGMENTED_UNSAT_WITHOUT_PROOF":
            raise RuntimeError(f"unexpected configuration step {step.result}")
        selected = set(step.selected_indices)
        selected_order, _ = augmented_edges_from_pool(
            len(points),
            base_edges,
            candidates,
            candidate_edges,
            selected,
        )
        selected_points = [
            *points,
            *[candidates[index].point for index in selected_order],
        ]
        exact_edges = complete_unit_edges(selected_points)
        encoding = ColoringSAT(len(selected_points), exact_edges, color_count)
        cnf_path = output_path.with_suffix(".cnf")
        proof_path = output_path.with_suffix(".drat")
        graph_output_path = output_path.with_name(
            output_path.stem + "-selected.graph.json"
        )
        encoding.write_dimacs_cnf(cnf_path)
        proof_result = encoding.solve_with_drat(proof_path)
        if proof_result is not None:
            raise RuntimeError("proof solver disagrees with augmented UNSAT")
        proof_verified, transcript = _verify_drat(
            Path(proof_checker_path),
            cnf_path,
            proof_path,
        )
        _write_graph(
            graph_output_path,
            "Potential six-chromatic second-closure graph",
            selected_points,
            exact_edges,
        )
        report["potential_lower_bound_certificate"] = {
            "graph_path": str(graph_output_path),
            "cnf_path": str(cnf_path),
            "proof_path": str(proof_path),
            "vertices": len(selected_points),
            "exact_complete_edges": len(exact_edges),
            "drat_verified": proof_verified,
            "drat_transcript": transcript,
        }
        report["lower_bound_6_proved"] = proof_verified
        termination = (
            "EXACT_5_UNSAT_DRAT_VERIFIED"
            if proof_verified
            else "EXACT_5_UNSAT_PROOF_VERIFICATION_FAILED"
        )
        break

    report["status"] = termination
    report["resume_state"] = {
        "next_iteration": len(report["iterations"]),
        "selected_candidate_indices": sorted(selected),
        "current_base_coloring": list(current_coloring),
    }
    report["final_selected_candidates"] = len(selected)
    report["elapsed_seconds"] = previous_elapsed + time.monotonic() - started
    write_checkpoint(output_path, report)
    return report


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--graph",
        type=Path,
        default=Path(
            "artifacts/cegis/"
            "heule529-degree3-configurations-full-pool.graph.json"
        ),
    )
    parser.add_argument(
        "--coloring-report",
        type=Path,
        default=Path(
            "artifacts/cegis/heule529-degree3-configurations.json"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/cegis/heule1061-second-closure.json"),
    )
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--min-degree", type=int, default=3)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--finalize-checkpoint", action="store_true")
    parser.add_argument("--direct-audit", action="store_true")
    parser.add_argument("--complete-audit", action="store_true")
    parser.add_argument(
        "--proof-checker",
        type=Path,
        default=Path("third_party/drat-trim/drat-trim"),
    )
    args = parser.parse_args()
    if args.direct_audit:
        report = run_direct_closure_audit(
            args.graph,
            args.output,
            minimum_candidate_degree=args.min_degree,
            proof_checker_path=args.proof_checker,
        )
    elif args.finalize_checkpoint:
        report = finalize_closure_checkpoint(args.graph, args.output)
    else:
        report = run_closure_cegis(
            args.graph,
            args.coloring_report,
            args.output,
            max_iterations=args.iterations,
            minimum_candidate_degree=args.min_degree,
            resume=args.resume,
            complete_audit=args.complete_audit,
            proof_checker_path=args.proof_checker,
        )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "status": report["status"],
                "iterations": len(report.get("iterations", [])),
                "candidate_generation": report["candidate_generation"],
                "candidate_interactions": report["candidate_interactions"],
                "final_selected_candidates": report.get(
                    "final_selected_candidates", 0
                ),
                "lower_bound_6_proved": report["lower_bound_6_proved"],
                "elapsed_seconds": report["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
