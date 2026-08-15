"""Seeded TabuCol heuristic for finding (never proving) graph colorings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import time

import numpy as np

from coloring_sat import ColoringSAT, parse_dimacs_edge


def tabu_repair(
    vertex_count: int,
    edges: set[tuple[int, int]],
    initial: list[int],
    *,
    colors: int = 5,
    seed: int = 20260804,
    iterations: int = 2_000_000,
    sample_size: int = 256,
) -> tuple[list[int] | None, dict[str, object]]:
    """Try to remove monochromatic edges; return only a validated coloring."""

    rng = random.Random(seed)
    adjacency: list[list[int]] = [[] for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    coloring = np.asarray(initial, dtype=np.int8)
    counts = np.zeros((vertex_count, colors), dtype=np.int16)
    for left, right in edges:
        counts[left, coloring[right]] += 1
        counts[right, coloring[left]] += 1
    conflict_degree = counts[np.arange(vertex_count), coloring].astype(np.int32)
    total = int(conflict_degree.sum() // 2)
    best_total = total
    best_coloring = coloring.copy()
    tabu_until = np.zeros((vertex_count, colors), dtype=np.int64)
    active = [vertex for vertex in range(vertex_count) if conflict_degree[vertex]]
    position = np.full(vertex_count, -1, dtype=np.int32)
    for index, vertex in enumerate(active):
        position[vertex] = index

    def refresh(vertex: int) -> None:
        is_active = bool(conflict_degree[vertex])
        index = int(position[vertex])
        if is_active and index < 0:
            position[vertex] = len(active)
            active.append(vertex)
        elif not is_active and index >= 0:
            last = active.pop()
            if index < len(active):
                active[index] = last
                position[last] = index
            position[vertex] = -1

    started = time.monotonic()
    for iteration in range(1, iterations + 1):
        if total == 0:
            result = [int(color) for color in coloring]
            if not ColoringSAT(
                vertex_count,
                edges,
                colors,
                break_color_symmetry=False,
            ).validate(result):
                raise AssertionError("TabuCol produced an invalid coloring")
            return result, {
                "status": "SAT",
                "seed": seed,
                "iterations": iteration - 1,
                "initial_conflicts": int(
                    sum(initial[left] == initial[right] for left, right in edges)
                ),
                "best_conflicts": 0,
                "elapsed_seconds": time.monotonic() - started,
            }
        pool = (
            active
            if len(active) <= sample_size
            else rng.sample(active, sample_size)
        )
        best_moves: list[tuple[int, int, int]] = []
        best_delta = 10**9
        for vertex in pool:
            old = int(coloring[vertex])
            old_conflicts = int(counts[vertex, old])
            for new in range(colors):
                if new == old:
                    continue
                delta = int(counts[vertex, new]) - old_conflicts
                forbidden = tabu_until[vertex, new] > iteration
                if forbidden and total + delta >= best_total:
                    continue
                if delta < best_delta:
                    best_delta = delta
                    best_moves = [(vertex, old, new)]
                elif delta == best_delta:
                    best_moves.append((vertex, old, new))
        if not best_moves:
            tabu_until.fill(0)
            continue
        vertex, old, new = rng.choice(best_moves)
        coloring[vertex] = new
        total += best_delta
        conflict_degree[vertex] = counts[vertex, new]
        refresh(vertex)
        for neighbor in adjacency[vertex]:
            counts[neighbor, old] -= 1
            counts[neighbor, new] += 1
            conflict_degree[neighbor] = counts[neighbor, coloring[neighbor]]
            refresh(neighbor)
        tenure = 7 + rng.randrange(11) + int(0.4 * len(active) ** 0.5)
        tabu_until[vertex, old] = iteration + tenure
        if total < best_total:
            best_total = total
            best_coloring = coloring.copy()
        if iteration % 250_000 == 0 and total > best_total + 100:
            coloring[:] = best_coloring
            for vertex_index in range(vertex_count):
                counts[vertex_index].fill(0)
            for left, right in edges:
                counts[left, coloring[right]] += 1
                counts[right, coloring[left]] += 1
            conflict_degree[:] = counts[np.arange(vertex_count), coloring]
            total = int(conflict_degree.sum() // 2)
            active[:] = [
                vertex for vertex in range(vertex_count) if conflict_degree[vertex]
            ]
            position.fill(-1)
            for index, active_vertex in enumerate(active):
                position[active_vertex] = index
            tabu_until.fill(0)
    return None, {
        "status": "UNKNOWN_ITERATION_LIMIT",
        "seed": seed,
        "iterations": iterations,
        "initial_conflicts": int(
            sum(initial[left] == initial[right] for left, right in edges)
        ),
        "best_conflicts": best_total,
        "elapsed_seconds": time.monotonic() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge", type=Path, required=True)
    parser.add_argument("--initial-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--iterations", type=int, default=2_000_000)
    parser.add_argument("--sample-size", type=int, default=256)
    args = parser.parse_args()
    vertex_count, edges = parse_dimacs_edge(args.edge)
    initial_payload = json.loads(args.initial_report.read_text(encoding="utf-8"))
    initial = [int(color) for color in initial_payload["whole_graph_coloring"]]
    coloring, stats = tabu_repair(
        vertex_count,
        edges,
        initial,
        seed=args.seed,
        iterations=args.iterations,
        sample_size=args.sample_size,
    )
    payload = {**stats, "coloring": coloring}
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({key: value for key, value in payload.items() if key != "coloring"}))
    return 10 if coloring is not None else 0


if __name__ == "__main__":
    raise SystemExit(main())
