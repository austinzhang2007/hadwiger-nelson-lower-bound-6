"""Exact SAT repair of a coloring inside a frozen graph neighbourhood.

An UNSAT result from this module is deliberately labelled local: it says only
that the chosen active neighbourhood cannot be recolored while every outside
vertex keeps its input color.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

from pysat.solvers import Solver

from coloring_sat import Edge, parse_dimacs_edge


def conflict_edges(
    edges: Iterable[Edge], coloring: Sequence[int]
) -> tuple[Edge, ...]:
    """Return normalized monochromatic edges in deterministic order."""

    normalized = {(min(left, right), max(left, right)) for left, right in edges}
    return tuple(
        edge
        for edge in sorted(normalized)
        if coloring[edge[0]] == coloring[edge[1]]
    )


def _active_neighbourhood(
    vertex_count: int,
    edges: set[Edge],
    seeds: set[int],
    radius: int,
) -> tuple[set[int], list[list[int]]]:
    adjacency: list[list[int]] = [[] for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    active = set(seeds)
    frontier = set(seeds)
    for _ in range(radius):
        following = {
            neighbor
            for vertex in frontier
            for neighbor in adjacency[vertex]
            if neighbor not in active
        }
        active.update(following)
        frontier = following
        if not frontier:
            break
    return active, adjacency


def repair_coloring_locally(
    vertex_count: int,
    edges: Iterable[Edge],
    initial: Sequence[int],
    *,
    colors: int = 5,
    radius: int = 0,
    solver_name: str = "cadical195",
) -> tuple[tuple[int, ...] | None, dict[str, int | str]]:
    """Recolor a conflict neighbourhood while freezing all other vertices."""

    if vertex_count < 0 or colors < 1 or radius < 0:
        raise ValueError("invalid vertex, color, or radius value")
    if len(initial) != vertex_count or any(
        color < 0 or color >= colors for color in initial
    ):
        raise ValueError("initial coloring has invalid length or color")
    normalized = {(min(left, right), max(left, right)) for left, right in edges}
    if any(
        left < 0 or right >= vertex_count or left == right
        for left, right in normalized
    ):
        raise ValueError("edge endpoint outside graph")
    initial_conflicts = conflict_edges(normalized, initial)
    if not initial_conflicts:
        return tuple(initial), {
            "status": "SAT",
            "radius": radius,
            "active_vertices": 0,
            "initial_conflicts": 0,
            "final_conflicts": 0,
        }

    seeds = {vertex for edge in initial_conflicts for vertex in edge}
    active, adjacency = _active_neighbourhood(
        vertex_count, normalized, seeds, radius
    )
    ordered = sorted(active)
    local_index = {vertex: index for index, vertex in enumerate(ordered)}

    def variable(vertex: int, color: int) -> int:
        return local_index[vertex] * colors + color + 1

    clauses: list[list[int]] = []
    for vertex in ordered:
        variables = [variable(vertex, color) for color in range(colors)]
        clauses.append(variables)
        for left_color in range(colors):
            for right_color in range(left_color + 1, colors):
                clauses.append(
                    [
                        -variable(vertex, left_color),
                        -variable(vertex, right_color),
                    ]
                )
        forbidden = {
            initial[neighbor]
            for neighbor in adjacency[vertex]
            if neighbor not in active
        }
        clauses.extend([[-variable(vertex, color)] for color in forbidden])
    for left, right in normalized:
        if left in active and right in active:
            clauses.extend(
                [
                    -variable(left, color),
                    -variable(right, color),
                ]
                for color in range(colors)
            )

    with Solver(name=solver_name, bootstrap_with=clauses) as solver:
        solver.set_phases(
            [variable(vertex, initial[vertex]) for vertex in ordered]
        )
        if not solver.solve():
            return None, {
                "status": "UNKNOWN_LOCAL_UNSAT",
                "radius": radius,
                "active_vertices": len(active),
                "initial_conflicts": len(initial_conflicts),
                "final_conflicts": len(initial_conflicts),
            }
        positive = {literal for literal in solver.get_model() if literal > 0}

    repaired = list(initial)
    for vertex in ordered:
        selected = [
            color
            for color in range(colors)
            if variable(vertex, color) in positive
        ]
        if len(selected) != 1:
            raise RuntimeError("local SAT model does not select one color")
        repaired[vertex] = selected[0]
    remaining = conflict_edges(normalized, repaired)
    if remaining:
        raise RuntimeError("local SAT model failed full-edge validation")
    return tuple(repaired), {
        "status": "SAT",
        "radius": radius,
        "active_vertices": len(active),
        "initial_conflicts": len(initial_conflicts),
        "final_conflicts": 0,
    }


def _load_coloring(path: Path) -> list[int]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for key in ("whole_graph_coloring", "coloring", "best_coloring"):
        raw = payload.get(key)
        if isinstance(raw, list):
            return [int(color) for color in raw]
    raise ValueError("initial report has no coloring")


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge", type=Path, required=True)
    parser.add_argument("--initial-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--colors", type=int, default=5)
    parser.add_argument("--radius", type=int, default=0)
    parser.add_argument("--solver", default="cadical195")
    args = parser.parse_args()

    vertex_count, edges = parse_dimacs_edge(args.edge)
    coloring, stats = repair_coloring_locally(
        vertex_count,
        edges,
        _load_coloring(args.initial_report),
        colors=args.colors,
        radius=args.radius,
        solver_name=args.solver,
    )
    payload: dict[str, object] = {**stats, "coloring": coloring}
    _write_json(args.output, payload)
    print(json.dumps(stats, sort_keys=True))
    return 10 if coloring is not None else 0


if __name__ == "__main__":
    raise SystemExit(main())
