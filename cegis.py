"""Coloring-blocking CEGIS for exact unit-distance graph augmentation."""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Hashable, Iterable, Mapping, Sequence, TypeVar

import numpy as np
import sympy as sp
from scipy.spatial import cKDTree
from pysat.solvers import Solver

from coloring_sat import ColoringSAT, parse_dimacs_edge
from exact_geometry import (
    Point,
    complete_bipartite_unit_edges,
    complete_unit_edges,
    parse_mathematica_vertices,
    serialize_graph,
)
from exact_geometry import unit_circle_intersections


T = TypeVar("T", bound=Hashable)


@dataclass(frozen=True)
class CandidateStats:
    center_pairs: int
    raw_intersections: int
    distinct_numeric_intersections: int
    numeric_degree_candidates: int
    certified_candidates: int
    elapsed_seconds: float


@dataclass(frozen=True)
class CertifiedCandidate:
    point: Point
    neighbors: tuple[int, ...]
    generator_pair: tuple[int, int]


@dataclass(frozen=True)
class _NumericSeed:
    x: float
    y: float
    generator_pair: tuple[int, int]
    possible_neighbors: tuple[int, ...]


def _numeric_key(x: float, y: float, digits: int) -> tuple[float, float]:
    return round(x, digits), round(y, digits)


def _numeric_candidate_seeds(
    points: Sequence[Point],
    *,
    min_neighbors: int,
    tolerance: float,
    dedup_digits: int,
) -> tuple[list[_NumericSeed], tuple[int, int, int]]:
    coordinates = np.asarray([point.approximate() for point in points], dtype=float)
    if len(coordinates) < 2:
        return [], (0, 0, 0)
    tree = cKDTree(coordinates)
    center_pairs = sorted(tree.query_pairs(r=2.0 + tolerance))
    existing_keys = {
        _numeric_key(float(x), float(y), dedup_digits) for x, y in coordinates
    }
    intersections: dict[
        tuple[float, float], tuple[float, float, tuple[int, int]]
    ] = {}
    raw_intersections = 0
    for i, j in center_pairs:
        dx, dy = coordinates[j] - coordinates[i]
        d2 = float(dx * dx + dy * dy)
        if d2 <= tolerance * tolerance or d2 > 4.0 + tolerance:
            continue
        distance = math.sqrt(max(d2, 0.0))
        height = math.sqrt(max(0.0, 1.0 - d2 / 4.0))
        midpoint = (coordinates[i] + coordinates[j]) / 2.0
        perpendicular = np.asarray([-dy / distance, dx / distance])
        signs = (1.0,) if height <= tolerance else (-1.0, 1.0)
        for sign in signs:
            raw_intersections += 1
            candidate = midpoint + sign * height * perpendicular
            x, y = float(candidate[0]), float(candidate[1])
            key = _numeric_key(x, y, dedup_digits)
            if key in existing_keys:
                continue
            intersections.setdefault(key, (x, y, (i, j)))

    seeds: list[_NumericSeed] = []
    for x, y, generator_pair in intersections.values():
        possible = tree.query_ball_point([x, y], r=1.0 + tolerance)
        exact_radius_numerically = tuple(
            sorted(
                index
                for index in possible
                if abs(
                    (coordinates[index, 0] - x) ** 2
                    + (coordinates[index, 1] - y) ** 2
                    - 1.0
                )
                <= 4.0 * tolerance
            )
        )
        if len(exact_radius_numerically) >= min_neighbors:
            seeds.append(
                _NumericSeed(x, y, generator_pair, exact_radius_numerically)
            )
    return seeds, (len(center_pairs), raw_intersections, len(intersections))


def _recover_seed(seed: _NumericSeed, points: Sequence[Point]) -> CertifiedCandidate:
    i, j = seed.generator_pair
    intersections = unit_circle_intersections(points[i], points[j])
    if not intersections:
        raise RuntimeError("numeric circle intersection has no exact counterpart")
    point = min(
        intersections,
        key=lambda candidate: (
            candidate.approximate()[0] - seed.x
        )
        ** 2
        + (candidate.approximate()[1] - seed.y) ** 2,
    )
    neighbors = tuple(
        index
        for index in seed.possible_neighbors
        if point.squared_distance(points[index]) == 1
    )
    return CertifiedCandidate(point, neighbors, seed.generator_pair)


def generate_certified_unit_circle_candidates(
    points: Sequence[Point],
    *,
    min_neighbors: int = 5,
    tolerance: float = 1e-9,
    dedup_digits: int = 14,
) -> tuple[list[CertifiedCandidate], CandidateStats]:
    """Numerically prefilter all unit-circle intersections, then certify exactly."""

    started = time.monotonic()
    seeds, counts = _numeric_candidate_seeds(
        points,
        min_neighbors=min_neighbors,
        tolerance=tolerance,
        dedup_digits=dedup_digits,
    )
    certified_by_point: dict[Point, CertifiedCandidate] = {}
    existing_points = set(points)
    for seed in seeds:
        candidate = _recover_seed(seed, points)
        if (
            candidate.point in existing_points
            or len(candidate.neighbors) < min_neighbors
        ):
            continue
        # Every recorded relation is exact; numeric prefiltering can only cause
        # a candidate to be omitted, never a false unit edge to be accepted.
        if not all(
            candidate.point.squared_distance(points[index]) == 1
            for index in candidate.neighbors
        ):
            raise RuntimeError("uncertified candidate neighbor escaped recovery")
        previous = certified_by_point.get(candidate.point)
        if previous is None:
            certified_by_point[candidate.point] = candidate
        else:
            certified_by_point[candidate.point] = CertifiedCandidate(
                point=candidate.point,
                neighbors=tuple(
                    sorted(set(previous.neighbors).union(candidate.neighbors))
                ),
                generator_pair=previous.generator_pair,
            )
    center_pairs, raw_intersections, distinct_intersections = counts
    candidates = sorted(
        certified_by_point.values(),
        key=lambda item: (
            -len(item.neighbors),
            item.point.approximate()[0],
            item.point.approximate()[1],
        ),
    )
    candidates = [
        candidate
        for candidate in recertify_candidate_neighbors(points, candidates)
        if len(candidate.neighbors) >= min_neighbors
    ]
    candidates.sort(
        key=lambda item: (
            -len(item.neighbors),
            item.point.approximate()[0],
            item.point.approximate()[1],
        )
    )
    stats = CandidateStats(
        center_pairs=center_pairs,
        raw_intersections=raw_intersections,
        distinct_numeric_intersections=distinct_intersections,
        numeric_degree_candidates=len(seeds),
        certified_candidates=len(candidates),
        elapsed_seconds=time.monotonic() - started,
    )
    return candidates, stats


def recertify_candidate_neighbors(
    base_points: Sequence[Point],
    candidates: Sequence[CertifiedCandidate],
) -> list[CertifiedCandidate]:
    """Replace numeric neighbor hints with the complete exact cross adjacency."""

    candidates = list(candidates)
    cross_edges = complete_bipartite_unit_edges(
        base_points,
        [candidate.point for candidate in candidates],
    )
    neighbors: list[list[int]] = [[] for _ in candidates]
    for base_index, candidate_index in sorted(cross_edges):
        neighbors[candidate_index].append(base_index)
    return [
        CertifiedCandidate(
            point=candidate.point,
            neighbors=tuple(neighbors[index]),
            generator_pair=candidate.generator_pair,
        )
        for index, candidate in enumerate(candidates)
    ]


def blocking_candidates(
    candidate_neighbors: Mapping[T, Sequence[int]],
    colorings: Sequence[Sequence[int]],
    color_count: int,
) -> dict[T, frozenset[int]]:
    coverage: dict[T, frozenset[int]] = {}
    all_colors = set(range(color_count))
    for candidate, neighbors in candidate_neighbors.items():
        blocked = {
            coloring_index
            for coloring_index, coloring in enumerate(colorings)
            if {coloring[index] for index in neighbors} == all_colors
        }
        coverage[candidate] = frozenset(blocked)
    return coverage


def interacting_pair_coverage(
    candidate_neighbors: Mapping[T, Sequence[int]],
    candidate_unit_edges: Iterable[tuple[T, T]],
    colorings: Sequence[Sequence[int]],
    color_count: int,
) -> dict[tuple[T, T], frozenset[int]]:
    """Coverage of adjacent candidate pairs under full two-vertex coloring.

    A pair gives genuinely new blocking power when each point individually has
    the same single available color and the unit edge between them forbids both
    from taking it.
    """

    all_colors = set(range(color_count))
    coverage: dict[tuple[T, T], frozenset[int]] = {}
    for left, right in candidate_unit_edges:
        blocked: set[int] = set()
        for coloring_index, coloring in enumerate(colorings):
            left_allowed = all_colors - {
                coloring[index] for index in candidate_neighbors[left]
            }
            right_allowed = all_colors - {
                coloring[index] for index in candidate_neighbors[right]
            }
            feasible = any(
                left_color != right_color
                for left_color in left_allowed
                for right_color in right_allowed
            )
            if not feasible:
                blocked.add(coloring_index)
        coverage[(left, right)] = frozenset(blocked)
    return coverage


def candidate_extension(
    candidate_neighbors: Mapping[T, Sequence[int]],
    candidate_unit_edges: Iterable[tuple[T, T]],
    base_coloring: Sequence[int],
    color_count: int,
    *,
    subset: Iterable[T] | None = None,
) -> dict[T, int] | None:
    """Color a candidate configuration subject to colors forbidden by the base."""

    candidates = (
        list(candidate_neighbors)
        if subset is None
        else list(subset)
    )
    candidate_set = set(candidates)
    position = {candidate: index for index, candidate in enumerate(candidates)}

    def variable(candidate: T, color: int) -> int:
        return position[candidate] * color_count + color + 1

    clauses: list[list[int]] = []
    all_colors = set(range(color_count))
    for candidate in candidates:
        allowed = all_colors - {
            base_coloring[index] for index in candidate_neighbors[candidate]
        }
        if not allowed:
            return None
        clauses.append([variable(candidate, color) for color in sorted(allowed)])
        for color in all_colors - allowed:
            clauses.append([-variable(candidate, color)])
        for left in range(color_count):
            for right in range(left + 1, color_count):
                clauses.append(
                    [
                        -variable(candidate, left),
                        -variable(candidate, right),
                    ]
                )
    for left, right in candidate_unit_edges:
        if left not in candidate_set or right not in candidate_set:
            continue
        for color in range(color_count):
            clauses.append(
                [-variable(left, color), -variable(right, color)]
            )
    with Solver(name="cadical195", bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None
        positive = {literal for literal in solver.get_model() if literal > 0}
    return {
        candidate: next(
            color
            for color in range(color_count)
            if variable(candidate, color) in positive
        )
        for candidate in candidates
    }


def candidate_pool_extension_conclusion(candidate_count: int) -> str:
    return (
        f"No subset of this {candidate_count}-candidate pool can block "
        "this base coloring."
    )


def minimize_blocking_configuration(
    candidate_neighbors: Mapping[T, Sequence[int]],
    candidate_unit_edges: Iterable[tuple[T, T]],
    base_coloring: Sequence[int],
    color_count: int,
) -> tuple[T, ...] | None:
    """Return a deletion-minimal unextendable candidate set, or ``None``."""

    edges = list(candidate_unit_edges)
    core = list(candidate_neighbors)
    if candidate_extension(
        candidate_neighbors, edges, base_coloring, color_count, subset=core
    ) is not None:
        return None
    for candidate in list(core):
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


def greedy_cover(
    coverage: Mapping[T, frozenset[int]], universe: Iterable[int]
) -> tuple[list[T], set[int]]:
    uncovered = set(universe)
    available = list(coverage)
    chosen: list[T] = []
    while uncovered:
        best = None
        best_gain: set[int] = set()
        for candidate in available:
            gain = uncovered.intersection(coverage[candidate])
            if len(gain) > len(best_gain):
                best, best_gain = candidate, gain
        if best is None:
            break
        chosen.append(best)
        uncovered.difference_update(best_gain)
        available.remove(best)
    return chosen, uncovered


def _best_pair(
    coverage: Mapping[CertifiedCandidate, frozenset[int]],
    *,
    pool_size: int = 200,
) -> tuple[tuple[CertifiedCandidate, CertifiedCandidate] | None, frozenset[int]]:
    ranked = sorted(
        coverage,
        key=lambda candidate: len(coverage[candidate]),
        reverse=True,
    )[:pool_size]
    best_pair = None
    best_union: frozenset[int] = frozenset()
    for left_index, left in enumerate(ranked):
        for right in ranked[left_index + 1 :]:
            union = coverage[left] | coverage[right]
            if len(union) > len(best_union):
                best_pair, best_union = (left, right), union
    return best_pair, best_union


def _candidate_json(
    candidate: CertifiedCandidate, blocked: frozenset[int]
) -> dict[str, object]:
    x, y = candidate.point.approximate()
    return {
        "x_exact": sp.sstr(candidate.point.x),
        "y_exact": sp.sstr(candidate.point.y),
        "x_approx": x,
        "y_approx": y,
        "degree_to_base": len(candidate.neighbors),
        "neighbors_zero_based": list(candidate.neighbors),
        "generator_pair_zero_based": list(candidate.generator_pair),
        "blocked_coloring_indices": sorted(blocked),
        "blocked_count": len(blocked),
    }


def _augment_edges(
    base_vertex_count: int,
    base_edges: set[tuple[int, int]],
    candidates: Sequence[CertifiedCandidate],
) -> set[tuple[int, int]]:
    edges = set(base_edges)
    for offset, candidate in enumerate(candidates):
        vertex = base_vertex_count + offset
        edges.update((neighbor, vertex) for neighbor in candidate.neighbors)
    for left_index, left in enumerate(candidates):
        for right_index in range(left_index + 1, len(candidates)):
            right = candidates[right_index]
            if left.point.squared_distance(right.point) == 1:
                edges.add(
                    (
                        base_vertex_count + left_index,
                        base_vertex_count + right_index,
                    )
                )
    return edges


def run_phase1_cegis(
    coordinate_path: str | Path,
    edge_path: str | Path,
    output_path: str | Path,
    *,
    coloring_limit: int = 128,
    seed: int = 20260725,
    max_iterations: int = 5,
    minimum_candidate_degree: int = 4,
) -> dict[str, object]:
    started = time.monotonic()
    points = parse_mathematica_vertices(coordinate_path)
    vertex_count, base_edges = parse_dimacs_edge(edge_path)
    if vertex_count != len(points):
        raise ValueError("coordinate/edge vertex counts disagree")

    colorings = list(
        ColoringSAT(vertex_count, base_edges, 5).enumerate_colorings(
            coloring_limit
        )
    )
    if not colorings:
        raise RuntimeError("base graph unexpectedly has no 5-coloring")
    initial_coloring_count = len(colorings)

    candidates, candidate_stats = generate_certified_unit_circle_candidates(
        points, min_neighbors=minimum_candidate_degree
    )
    neighbors = {index: candidate.neighbors for index, candidate in enumerate(candidates)}
    candidate_unit_edges = sorted(
        complete_unit_edges([candidate.point for candidate in candidates])
    )
    initial_coverage = blocking_candidates(neighbors, colorings, 5)
    if not initial_coverage or not max(map(len, initial_coverage.values())):
        raise RuntimeError("no unit-circle candidate blocks any sampled coloring")

    selected_indices: set[int] = set()
    iteration_reports: list[dict[str, object]] = []
    stagnation_analysis: list[dict[str, object]] = []
    termination = "ITERATION_LIMIT"
    for iteration in range(max_iterations):
        single_coverage = blocking_candidates(neighbors, colorings, 5)
        pair_coverage = interacting_pair_coverage(
            neighbors, candidate_unit_edges, colorings, 5
        )
        chosen_single_indices, single_uncovered = greedy_cover(
            single_coverage, universe=range(len(colorings))
        )
        chosen_blockers: list[tuple[int, ...]] = [
            (index,) for index in chosen_single_indices
        ]
        uncovered = set(single_uncovered)
        if uncovered:
            chosen_pairs, uncovered = greedy_cover(
                pair_coverage, universe=uncovered
            )
            chosen_blockers.extend(chosen_pairs)
        if uncovered:
            termination = "SINGLE_AND_INTERACTING_PAIR_STAGNATION"
            for coloring_index in sorted(uncovered):
                extension = candidate_extension(
                    neighbors,
                    candidate_unit_edges,
                    colorings[coloring_index],
                    5,
                )
                if extension is not None:
                    full_coloring = [
                        *colorings[coloring_index],
                        *[
                            extension[index]
                            for index in range(len(candidates))
                        ],
                    ]
                    stagnation_analysis.append(
                        {
                            "coloring_index": coloring_index,
                            "entire_candidate_pool_extends": True,
                            "candidate_assignment": [
                                extension[index]
                                for index in range(len(candidates))
                            ],
                            "full_pool_coloring": full_coloring,
                            "conclusion": candidate_pool_extension_conclusion(
                                len(candidates)
                            ),
                        }
                    )
                else:
                    core = minimize_blocking_configuration(
                        neighbors,
                        candidate_unit_edges,
                        colorings[coloring_index],
                        5,
                    )
                    if core is None:
                        raise RuntimeError("UNSAT extension lost during core extraction")
                    stagnation_analysis.append(
                        {
                            "coloring_index": coloring_index,
                            "entire_candidate_pool_extends": False,
                            "minimal_blocking_candidate_indices": list(core),
                            "minimal_blocking_size": len(core),
                            "minimal_blocking_candidates": [
                                _candidate_json(
                                    candidates[index],
                                    single_coverage[index],
                                )
                                for index in core
                            ],
                        }
                    )
            iteration_reports.append(
                {
                    "iteration": iteration,
                    "known_colorings": len(colorings),
                    "result": termination,
                    "uncovered_coloring_indices": sorted(uncovered),
                }
            )
            break
        additions = {
            candidate_index
            for blocker in chosen_blockers
            for candidate_index in blocker
        } - selected_indices
        if not additions:
            raise RuntimeError("CEGIS cover added no point for a new coloring")
        selected_indices.update(additions)
        selected_candidates = [
            candidates[index] for index in sorted(selected_indices)
        ]
        augmented_edges = _augment_edges(
            vertex_count, base_edges, selected_candidates
        )
        new_coloring = ColoringSAT(
            vertex_count + len(selected_candidates), augmented_edges, 5
        ).solve()
        iteration_report: dict[str, object] = {
            "iteration": iteration,
            "known_colorings_before_solve": len(colorings),
            "chosen_blockers": [list(blocker) for blocker in chosen_blockers],
            "new_candidate_indices": sorted(additions),
            "selected_candidate_indices": sorted(selected_indices),
            "augmented_vertices": vertex_count + len(selected_candidates),
            "augmented_edges": len(augmented_edges),
        }
        if new_coloring is None:
            termination = "UNSAT_WITHOUT_PROOF_NOT_A_CERTIFICATE"
            iteration_report["result"] = termination
            iteration_reports.append(iteration_report)
            break
        base_coloring = tuple(new_coloring[:vertex_count])
        if base_coloring in set(colorings):
            raise RuntimeError("CEGIS returned an already-blocked base coloring")
        colorings.append(base_coloring)
        candidate_color_details = []
        for position, candidate_index in enumerate(sorted(selected_indices)):
            candidate = candidates[candidate_index]
            neighbor_colors = sorted(
                {new_coloring[index] for index in candidate.neighbors}
            )
            assigned_color = new_coloring[vertex_count + position]
            if assigned_color in neighbor_colors:
                raise RuntimeError("SAT model violates an augmented unit edge")
            candidate_color_details.append(
                {
                    "candidate_index": candidate_index,
                    "assigned_color": assigned_color,
                    "neighbor_colors": neighbor_colors,
                }
            )
        iteration_report.update(
            {
                "result": "SAT_COUNTEREXAMPLE",
                "new_coloring_index": len(colorings) - 1,
                "candidate_color_details": candidate_color_details,
                "new_base_coloring": list(base_coloring),
            }
        )
        iteration_reports.append(iteration_report)

    final_single_coverage = blocking_candidates(neighbors, colorings, 5)
    final_pair_coverage = interacting_pair_coverage(
        neighbors, candidate_unit_edges, colorings, 5
    )
    ranked_indices = sorted(
        range(len(candidates)),
        key=lambda index: (
            -len(final_single_coverage[index]),
            -len(candidates[index].neighbors),
            sp.sstr(candidates[index].point.x),
            sp.sstr(candidates[index].point.y),
        ),
    )
    final_single_chosen, final_single_uncovered = greedy_cover(
        final_single_coverage, universe=range(len(colorings))
    )
    final_pair_chosen: list[tuple[int, int]] = []
    final_uncovered = set(final_single_uncovered)
    if final_uncovered:
        final_pair_chosen, final_uncovered = greedy_cover(
            final_pair_coverage, universe=final_uncovered
        )
    pair_joint_only = {
        pair: blocked
        - (
            final_single_coverage[pair[0]]
            | final_single_coverage[pair[1]]
        )
        for pair, blocked in final_pair_coverage.items()
    }
    best_pair_indices = max(
        final_pair_coverage,
        key=lambda pair: (
            len(pair_joint_only[pair]),
            len(final_pair_coverage[pair]),
        ),
        default=None,
    )
    best_pair_json: dict[str, object] | None = None
    if best_pair_indices is not None:
        left, right = best_pair_indices
        joint_only = pair_joint_only[best_pair_indices]
        best_pair_json = {
            "candidate_indices": [left, right],
            "blocked_count": len(final_pair_coverage[best_pair_indices]),
            "joint_only_blocked_indices": sorted(joint_only),
            "joint_only_blocked_count": len(joint_only),
            "candidates": [
                _candidate_json(candidates[index], final_single_coverage[index])
                for index in best_pair_indices
            ],
        }

    report: dict[str, object] = {
        "status": "NUMERIC_PREFILTER_EXACT_RETAINED_CANDIDATES",
        "lower_bound_6_proved": False,
        "seed": seed,
        "base": {
            "vertices": vertex_count,
            "edges": len(base_edges),
            "coordinate_source": str(coordinate_path),
            "edge_source": str(edge_path),
        },
        "initial_canonical_5_colorings": initial_coloring_count,
        "total_known_5_colorings_after_cegis": len(colorings),
        "candidate_generation": asdict(candidate_stats),
        "minimum_candidate_degree": minimum_candidate_degree,
        "candidate_pool": [
            _candidate_json(
                candidate, final_single_coverage[index]
            )
            for index, candidate in enumerate(candidates)
        ],
        "candidate_unit_edges_zero_based": [
            list(edge) for edge in candidate_unit_edges
        ],
        "blocking": {
            "candidates_blocking_at_least_one": sum(
                bool(blocked) for blocked in final_single_coverage.values()
            ),
            "candidate_interaction_edges": len(candidate_unit_edges),
            "best": _candidate_json(
                candidates[ranked_indices[0]],
                final_single_coverage[ranked_indices[0]],
            ),
            "top_20": [
                _candidate_json(
                    candidates[index], final_single_coverage[index]
                )
                for index in ranked_indices[:20]
            ],
            "greedy_single_candidate_indices": final_single_chosen,
            "greedy_uncovered_by_singles": sorted(final_single_uncovered),
            "greedy_interacting_pairs_for_single_uncovered": [
                list(pair) for pair in final_pair_chosen
            ],
            "greedy_uncovered_indices": sorted(final_uncovered),
            "interacting_pairs_with_joint_only_coverage": sum(
                bool(blocked) for blocked in pair_joint_only.values()
            ),
            "max_joint_only_pair_coverage": max(
                map(len, pair_joint_only.values()), default=0
            ),
            "best_interacting_pair": best_pair_json,
        },
        "cegis": {
            "max_iterations": max_iterations,
            "termination": termination,
            "selected_candidate_indices": sorted(selected_indices),
            "stagnation_analysis": stagnation_analysis,
            "iterations": iteration_reports,
        },
        "elapsed_seconds": time.monotonic() - started,
    }
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    full_pool_path = output.with_name(
        output.stem + "-full-pool.graph.json"
    )
    full_pool_edges = _augment_edges(vertex_count, base_edges, candidates)
    with full_pool_path.open("w", encoding="utf-8") as handle:
        json.dump(
            serialize_graph(
                "Heule 529 plus all certified unit-circle candidates",
                [*points, *[candidate.point for candidate in candidates]],
                full_pool_edges,
            ),
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
    report["full_candidate_pool_graph"] = {
        "path": str(full_pool_path),
        "vertices": vertex_count + len(candidates),
        "edges": len(full_pool_edges),
    }
    with output.open("w", encoding="utf-8") as handle:
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
        "--output",
        type=Path,
        default=Path("artifacts/cegis/heule529-unit-circles.json"),
    )
    parser.add_argument("--colorings", type=int, default=128)
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--min-degree", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260725)
    args = parser.parse_args()
    report = run_phase1_cegis(
        args.coordinates,
        args.edge_file,
        args.output,
        coloring_limit=args.colorings,
        seed=args.seed,
        max_iterations=args.iterations,
        minimum_candidate_degree=args.min_degree,
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "initial_colorings": report["initial_canonical_5_colorings"],
                "known_colorings": report["total_known_5_colorings_after_cegis"],
                "candidate_generation": report["candidate_generation"],
                "best": report["blocking"]["best"],
                "cegis_termination": report["cegis"]["termination"],
                "selected_candidate_indices": report["cegis"][
                    "selected_candidate_indices"
                ],
                "lower_bound_6_proved": report["lower_bound_6_proved"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
