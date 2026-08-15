"""CEGIS over interacting exact unit-circle candidate configurations."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Hashable, Iterable, Mapping, Sequence, TypeVar

import numpy as np
import sympy as sp
from scipy.spatial import cKDTree

from cegis import (
    CertifiedCandidate,
    candidate_extension,
    generate_certified_unit_circle_candidates,
)
from coloring_sat import ColoringSAT, parse_dimacs_edge
from exact_geometry import (
    Point,
    complete_unit_edges,
    parse_mathematica_vertices,
    parse_real_expr,
    serialize_graph,
)


T = TypeVar("T", bound=Hashable)


@dataclass(frozen=True)
class CandidateEdgeStats:
    numeric_pairs: int
    numeric_near_unit_pairs: int
    exact_unit_edges: int


@dataclass(frozen=True)
class ConfigurationStep:
    result: str
    blocking_core: tuple[int, ...]
    selected_indices: tuple[int, ...]
    augmented_edges: int
    new_base_coloring: tuple[int, ...] | None
    full_pool_assignment: tuple[int, ...] | None


def certified_candidate_unit_edges(
    candidates: Sequence[CertifiedCandidate],
    *,
    tolerance: float = 1e-9,
) -> tuple[set[tuple[int, int]], CandidateEdgeStats]:
    """Prefilter candidate pairs numerically and certify retained edges exactly."""

    if len(candidates) < 2:
        return set(), CandidateEdgeStats(0, 0, 0)
    coordinates = np.asarray(
        [candidate.point.approximate() for candidate in candidates],
        dtype=float,
    )
    possible = sorted(cKDTree(coordinates).query_pairs(r=1.0 + tolerance))
    near_unit = [
        (left, right)
        for left, right in possible
        if abs(
            float(np.sum((coordinates[left] - coordinates[right]) ** 2)) - 1.0
        )
        <= 4.0 * tolerance
    ]
    edges = {
        (left, right)
        for left, right in near_unit
        if candidates[left].point.squared_distance(candidates[right].point) == 1
    }
    return edges, CandidateEdgeStats(
        numeric_pairs=len(possible),
        numeric_near_unit_pairs=len(near_unit),
        exact_unit_edges=len(edges),
    )


def deletion_minimal_blocker(
    candidate_neighbors: Mapping[T, Sequence[int]],
    candidate_edges: Iterable[tuple[T, T]],
    base_coloring: Sequence[int],
    color_count: int,
) -> tuple[T, ...] | None:
    """Return a deletion-minimal unextendable candidate configuration."""

    edges = tuple(candidate_edges)
    core = list(candidate_neighbors)
    if candidate_extension(
        candidate_neighbors,
        edges,
        base_coloring,
        color_count,
        subset=core,
    ) is not None:
        return None
    for candidate in tuple(core):
        trial = [item for item in core if item != candidate]
        if candidate_extension(
            candidate_neighbors,
            edges,
            base_coloring,
            color_count,
            subset=trial,
        ) is None:
            core = trial
    return tuple(core)


def map_selected_candidates(
    selected_indices: Iterable[int],
    old_candidates: Sequence[CertifiedCandidate],
    expanded_candidates: Sequence[CertifiedCandidate],
) -> set[int]:
    """Map saved candidate indices to a new pool by exact coordinates."""

    expanded_by_point = {
        candidate.point: index for index, candidate in enumerate(expanded_candidates)
    }
    mapped: set[int] = set()
    for old_index in selected_indices:
        point = old_candidates[old_index].point
        if point not in expanded_by_point:
            raise ValueError(
                f"selected candidate {old_index} is absent from expanded candidate pool"
            )
        mapped.add(expanded_by_point[point])
    return mapped


def _candidate_from_report(record: Mapping[str, object]) -> CertifiedCandidate:
    return CertifiedCandidate(
        point=Point(
            parse_real_expr(str(record["x_exact"])),
            parse_real_expr(str(record["y_exact"])),
        ),
        neighbors=tuple(int(index) for index in record["neighbors_zero_based"]),
        generator_pair=tuple(
            int(index) for index in record["generator_pair_zero_based"]
        ),
    )


def resume_phase1_state(
    report: Mapping[str, object],
    expanded_candidates: Sequence[CertifiedCandidate],
    *,
    base_vertex_count: int,
) -> tuple[set[int], tuple[int, ...]]:
    """Recover selected points and the last stagnant base coloring exactly."""

    old_candidates = [
        _candidate_from_report(record) for record in report["candidate_pool"]
    ]
    cegis_record = report["cegis"]
    selected = map_selected_candidates(
        cegis_record["selected_candidate_indices"],
        old_candidates,
        expanded_candidates,
    )
    stagnation = cegis_record["stagnation_analysis"]
    if not stagnation:
        raise ValueError("phase-1 report has no stagnation coloring")
    full_coloring = stagnation[-1]["full_pool_coloring"]
    if len(full_coloring) < base_vertex_count:
        raise ValueError("phase-1 stagnation coloring is shorter than the base graph")
    return selected, tuple(
        int(color) for color in full_coloring[:base_vertex_count]
    )


def augmented_edges_from_pool(
    base_vertex_count: int,
    base_edges: Iterable[tuple[int, int]],
    candidates: Sequence[CertifiedCandidate],
    candidate_edges: Iterable[tuple[int, int]],
    selected_indices: Iterable[int],
) -> tuple[tuple[int, ...], set[tuple[int, int]]]:
    """Build an augmented graph using sorted candidate-pool indices."""

    order = tuple(sorted(set(selected_indices)))
    position = {
        pool_index: base_vertex_count + offset
        for offset, pool_index in enumerate(order)
    }
    edges = set(base_edges)
    for pool_index in order:
        vertex = position[pool_index]
        edges.update(
            (min(neighbor, vertex), max(neighbor, vertex))
            for neighbor in candidates[pool_index].neighbors
        )
    for left, right in candidate_edges:
        if left in position and right in position:
            u, v = position[left], position[right]
            edges.add((min(u, v), max(u, v)))
    return order, edges


def configuration_cegis_step(
    *,
    base_vertex_count: int,
    base_edges: Iterable[tuple[int, int]],
    candidates: Sequence[CertifiedCandidate],
    candidate_edges: Iterable[tuple[int, int]],
    selected_indices: Iterable[int],
    base_coloring: Sequence[int],
    color_count: int,
) -> ConfigurationStep:
    """Block one fixed base coloring or return its whole-pool extension."""

    candidate_edge_set = set(candidate_edges)
    neighbors = {
        index: candidate.neighbors for index, candidate in enumerate(candidates)
    }
    extension = candidate_extension(
        neighbors,
        candidate_edge_set,
        base_coloring,
        color_count,
    )
    selected = set(selected_indices)
    if extension is not None:
        return ConfigurationStep(
            result="WHOLE_POOL_EXTENDS",
            blocking_core=(),
            selected_indices=tuple(sorted(selected)),
            augmented_edges=0,
            new_base_coloring=None,
            full_pool_assignment=tuple(
                extension[index] for index in range(len(candidates))
            ),
        )

    core = deletion_minimal_blocker(
        neighbors,
        candidate_edge_set,
        base_coloring,
        color_count,
    )
    if core is None:
        raise RuntimeError("unextendable pool lost its blocking core")
    selected.update(core)
    order, augmented_edges = augmented_edges_from_pool(
        base_vertex_count,
        base_edges,
        candidates,
        candidate_edge_set,
        selected,
    )
    coloring = ColoringSAT(
        base_vertex_count + len(order),
        augmented_edges,
        color_count,
    ).solve()
    if coloring is None:
        return ConfigurationStep(
            result="AUGMENTED_UNSAT_WITHOUT_PROOF",
            blocking_core=tuple(core),
            selected_indices=order,
            augmented_edges=len(augmented_edges),
            new_base_coloring=None,
            full_pool_assignment=None,
        )
    return ConfigurationStep(
        result="SAT_COUNTEREXAMPLE",
        blocking_core=tuple(core),
        selected_indices=order,
        augmented_edges=len(augmented_edges),
        new_base_coloring=tuple(coloring[:base_vertex_count]),
        full_pool_assignment=None,
    )


def _candidate_record(
    candidate: CertifiedCandidate,
    base_coloring: Sequence[int] | None = None,
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
        neighbor_colors = sorted(
            {int(base_coloring[index]) for index in candidate.neighbors}
        )
        record["neighbor_colors"] = neighbor_colors
        record["allowed_colors"] = sorted(
            set(range(5)).difference(neighbor_colors)
        )
    return record


def _write_graph(
    path: Path,
    name: str,
    points: Sequence[Point],
    edges: Iterable[tuple[int, int]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            serialize_graph(name, points, edges),
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")


def _verify_drat(
    checker_path: Path,
    cnf_path: Path,
    proof_path: Path,
) -> tuple[bool, str]:
    completed = subprocess.run(
        [str(checker_path), str(cnf_path), str(proof_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    transcript = completed.stdout + completed.stderr
    return completed.returncode == 0 and "s VERIFIED" in transcript, transcript


def run_configuration_cegis(
    coordinate_path: str | Path,
    edge_path: str | Path,
    phase1_report_path: str | Path,
    output_path: str | Path,
    *,
    max_iterations: int = 100,
    minimum_candidate_degree: int = 3,
    color_count: int = 5,
    proof_checker_path: str | Path = "third_party/drat-trim/drat-trim",
) -> dict[str, object]:
    """Resume phase 1 with arbitrary interacting candidate configurations."""

    started = time.monotonic()
    coordinate_path = Path(coordinate_path)
    edge_path = Path(edge_path)
    phase1_report_path = Path(phase1_report_path)
    output_path = Path(output_path)
    points = parse_mathematica_vertices(coordinate_path)
    base_vertex_count, base_edges = parse_dimacs_edge(edge_path)
    if base_vertex_count != len(points):
        raise ValueError("coordinate/edge vertex counts disagree")
    candidates, candidate_stats = generate_certified_unit_circle_candidates(
        points,
        min_neighbors=minimum_candidate_degree,
    )
    candidate_edges, candidate_edge_stats = certified_candidate_unit_edges(
        candidates
    )
    with phase1_report_path.open(encoding="utf-8") as handle:
        phase1_report = json.load(handle)
    selected, current_base_coloring = resume_phase1_state(
        phase1_report,
        candidates,
        base_vertex_count=base_vertex_count,
    )
    base_encoding = ColoringSAT(
        base_vertex_count,
        base_edges,
        color_count,
    )
    if not base_encoding.validate(current_base_coloring):
        raise ValueError("phase-1 stagnation coloring is invalid on the base graph")

    report: dict[str, object] = {
        "status": "RUNNING",
        "lower_bound_6_proved": False,
        "base": {
            "vertices": base_vertex_count,
            "edges": len(base_edges),
            "coordinate_source": str(coordinate_path),
            "edge_source": str(edge_path),
        },
        "phase1_report": str(phase1_report_path),
        "minimum_candidate_degree": minimum_candidate_degree,
        "candidate_generation": asdict(candidate_stats),
        "candidate_interactions": asdict(candidate_edge_stats),
        "candidate_pool": [
            _candidate_record(candidate) for candidate in candidates
        ],
        "candidate_unit_edges_zero_based": [
            list(edge) for edge in sorted(candidate_edges)
        ],
        "initial_selected_candidate_indices": sorted(selected),
        "initial_selected_candidates": len(selected),
        "initial_base_coloring": list(current_base_coloring),
        "iterations": [],
    }
    termination = "ITERATION_LIMIT"
    strict_geometry_audit = False
    strict_coloring_audit = False
    for iteration_index in range(max_iterations):
        step_started = time.monotonic()
        fixed_coloring = current_base_coloring
        step = configuration_cegis_step(
            base_vertex_count=base_vertex_count,
            base_edges=base_edges,
            candidates=candidates,
            candidate_edges=candidate_edges,
            selected_indices=selected,
            base_coloring=fixed_coloring,
            color_count=color_count,
        )
        iteration: dict[str, object] = {
            "iteration": iteration_index,
            "result": step.result,
            "fixed_base_coloring": list(fixed_coloring),
            "blocking_core_candidate_indices": list(step.blocking_core),
            "blocking_core": [
                _candidate_record(candidates[index], fixed_coloring)
                for index in step.blocking_core
            ],
            "blocking_core_unit_edges": [
                list(edge)
                for edge in sorted(candidate_edges)
                if edge[0] in step.blocking_core
                and edge[1] in step.blocking_core
            ],
            "selected_candidate_indices": list(step.selected_indices),
            "selected_candidates": len(step.selected_indices),
            "augmented_edges": step.augmented_edges,
            "elapsed_seconds": time.monotonic() - step_started,
        }
        report["iterations"].append(iteration)

        if step.result == "WHOLE_POOL_EXTENDS":
            if step.full_pool_assignment is None:
                raise RuntimeError("whole-pool extension omitted its assignment")
            full_points = [
                *points,
                *[candidate.point for candidate in candidates],
            ]
            if len(set(full_points)) != len(full_points):
                raise RuntimeError("full candidate pool contains duplicate points")
            exact_full_edges = complete_unit_edges(full_points)
            strict_geometry_audit = True
            full_coloring = (
                *fixed_coloring,
                *step.full_pool_assignment,
            )
            strict_coloring_audit = ColoringSAT(
                len(full_points),
                exact_full_edges,
                color_count,
                break_color_symmetry=False,
            ).validate(full_coloring)
            full_pool_graph_path = output_path.with_name(
                output_path.stem + "-full-pool.graph.json"
            )
            _write_graph(
                full_pool_graph_path,
                "Heule 529 plus all degree-at-least-three candidates",
                full_points,
                exact_full_edges,
            )
            termination = (
                "EXACT_WHOLE_POOL_5_COLORING"
                if strict_coloring_audit
                else "PREFILTER_EXTENSION_FAILED_EXACT_EDGE_AUDIT"
            )
            report["whole_pool_failure_certificate"] = {
                "graph_path": str(full_pool_graph_path),
                "vertices": len(full_points),
                "exact_complete_edges": len(exact_full_edges),
                "coloring": list(full_coloring),
                "coloring_valid_on_exact_complete_edges": strict_coloring_audit,
            }
            break

        selected = set(step.selected_indices)
        if step.result == "SAT_COUNTEREXAMPLE":
            if step.new_base_coloring is None:
                raise RuntimeError("SAT counterexample omitted its base coloring")
            if step.new_base_coloring == fixed_coloring:
                raise RuntimeError("blocking configuration did not change the model")
            current_base_coloring = step.new_base_coloring
            iteration["new_base_coloring"] = list(current_base_coloring)
            continue

        if step.result != "AUGMENTED_UNSAT_WITHOUT_PROOF":
            raise RuntimeError(f"unexpected configuration step: {step.result}")
        selected_order, _ = augmented_edges_from_pool(
            base_vertex_count,
            base_edges,
            candidates,
            candidate_edges,
            selected,
        )
        selected_points = [
            *points,
            *[candidates[index].point for index in selected_order],
        ]
        if len(set(selected_points)) != len(selected_points):
            raise RuntimeError("selected graph contains duplicate points")
        exact_selected_edges = complete_unit_edges(selected_points)
        strict_geometry_audit = True
        encoding = ColoringSAT(
            len(selected_points),
            exact_selected_edges,
            color_count,
        )
        cnf_path = output_path.with_suffix(".cnf")
        proof_path = output_path.with_suffix(".drat")
        graph_path = output_path.with_name(
            output_path.stem + "-selected.graph.json"
        )
        encoding.write_dimacs_cnf(cnf_path)
        proof_result = encoding.solve_with_drat(proof_path)
        if proof_result is not None:
            raise RuntimeError("proof solver disagrees with augmented UNSAT result")
        proof_verified, proof_transcript = _verify_drat(
            Path(proof_checker_path),
            cnf_path,
            proof_path,
        )
        _write_graph(
            graph_path,
            "Potential six-chromatic exact unit-distance graph",
            selected_points,
            exact_selected_edges,
        )
        termination = (
            "EXACT_5_UNSAT_DRAT_VERIFIED"
            if proof_verified
            else "EXACT_5_UNSAT_PROOF_VERIFICATION_FAILED"
        )
        report["potential_lower_bound_certificate"] = {
            "graph_path": str(graph_path),
            "cnf_path": str(cnf_path),
            "proof_path": str(proof_path),
            "vertices": len(selected_points),
            "exact_complete_edges": len(exact_selected_edges),
            "drat_verified": proof_verified,
            "drat_transcript": proof_transcript,
        }
        report["lower_bound_6_proved"] = proof_verified
        break

    report["status"] = termination
    report["strict_geometry_audit"] = strict_geometry_audit
    report["strict_coloring_audit"] = strict_coloring_audit
    report["final_selected_candidate_indices"] = sorted(selected)
    report["final_selected_candidates"] = len(selected)
    report["elapsed_seconds"] = time.monotonic() - started
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coordinates",
        type=Path,
        default=Path("third_party/CNP-SAT/vtx/529.vtx"),
    )
    parser.add_argument(
        "--edge-file",
        type=Path,
        default=Path("third_party/CNP-SAT/edge/529.edge"),
    )
    parser.add_argument(
        "--phase1-report",
        type=Path,
        default=Path("artifacts/cegis/heule529-unit-circles.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/cegis/heule529-degree3-configurations.json"
        ),
    )
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--min-degree", type=int, default=3)
    parser.add_argument(
        "--proof-checker",
        type=Path,
        default=Path("third_party/drat-trim/drat-trim"),
    )
    args = parser.parse_args()
    report = run_configuration_cegis(
        args.coordinates,
        args.edge_file,
        args.phase1_report,
        args.output,
        max_iterations=args.iterations,
        minimum_candidate_degree=args.min_degree,
        proof_checker_path=args.proof_checker,
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "status": report["status"],
                "candidate_generation": report["candidate_generation"],
                "candidate_interactions": report["candidate_interactions"],
                "iterations": len(report["iterations"]),
                "final_selected_candidates": report[
                    "final_selected_candidates"
                ],
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
