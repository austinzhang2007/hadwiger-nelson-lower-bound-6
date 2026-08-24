"""Use SAT assumption cores to unfreeze only necessary coloring vertices.

This is a coloring search aid, not an UNSAT proof generator.  If the complete
formula becomes UNSAT, the status remains explicitly uncertified until a
separate DRAT/LRAT-producing run is checked.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

from pysat.solvers import Solver

from coloring_sat import ColoringSAT, Edge, parse_dimacs_edge
from local_coloring_repair import conflict_edges


def repair_coloring_core_guided(
    vertex_count: int,
    edges: Iterable[Edge],
    initial: Sequence[int],
    *,
    colors: int = 5,
    solver_name: str = "cadical195",
    max_rounds: int | None = None,
) -> tuple[tuple[int, ...] | None, dict[str, object]]:
    """Relax frozen colors in UNSAT cores until SAT or an uncertified stop."""

    if max_rounds is not None and max_rounds < 1:
        raise ValueError("max_rounds must be positive")
    normalized = {(min(left, right), max(left, right)) for left, right in edges}
    encoding = ColoringSAT(
        vertex_count,
        normalized,
        colors,
        break_color_symmetry=False,
        solver_name=solver_name,
    )
    if len(initial) != vertex_count or any(
        color < 0 or color >= colors for color in initial
    ):
        raise ValueError("initial coloring has invalid length or color")
    initial_bad = conflict_edges(normalized, initial)
    if not initial_bad:
        return tuple(initial), {
            "status": "SAT",
            "rounds": 0,
            "initial_conflicts": 0,
            "final_conflicts": 0,
            "relaxed_vertices": 0,
            "core_sizes": [],
        }

    relaxed = {vertex for edge in initial_bad for vertex in edge}
    core_sizes: list[int] = []
    rounds = 0
    phases = [
        encoding.variable(vertex, initial[vertex])
        for vertex in range(vertex_count)
    ]
    with Solver(
        name=solver_name, bootstrap_with=encoding.clauses
    ) as solver:
        solver.set_phases(phases)
        while max_rounds is None or rounds < max_rounds:
            rounds += 1
            assumptions = [
                encoding.variable(vertex, initial[vertex])
                for vertex in range(vertex_count)
                if vertex not in relaxed
            ]
            if solver.solve(assumptions=assumptions):
                model = {literal for literal in solver.get_model() if literal > 0}
                coloring = tuple(
                    next(
                        color
                        for color in range(colors)
                        if encoding.variable(vertex, color) in model
                    )
                    for vertex in range(vertex_count)
                )
                if not encoding.validate(coloring):
                    raise RuntimeError("core-guided solver returned invalid coloring")
                return coloring, {
                    "status": "SAT",
                    "rounds": rounds,
                    "initial_conflicts": len(initial_bad),
                    "final_conflicts": 0,
                    "relaxed_vertices": len(relaxed),
                    "core_sizes": core_sizes,
                }

            core = solver.get_core() or []
            core_vertices = {
                (abs(literal) - 1) // colors for literal in core
            }
            newly_relaxed = core_vertices - relaxed
            core_sizes.append(len(core_vertices))
            if not newly_relaxed:
                return None, {
                    "status": "UNKNOWN_FULL_UNSAT_UNCERTIFIED",
                    "rounds": rounds,
                    "initial_conflicts": len(initial_bad),
                    "final_conflicts": len(initial_bad),
                    "relaxed_vertices": len(relaxed),
                    "core_sizes": core_sizes,
                }
            relaxed.update(newly_relaxed)

    return None, {
        "status": "UNKNOWN_ROUND_LIMIT",
        "rounds": rounds,
        "initial_conflicts": len(initial_bad),
        "final_conflicts": len(initial_bad),
        "relaxed_vertices": len(relaxed),
        "core_sizes": core_sizes,
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
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--max-rounds", type=int)
    args = parser.parse_args()

    vertex_count, edges = parse_dimacs_edge(args.edge)
    coloring, stats = repair_coloring_core_guided(
        vertex_count,
        edges,
        _load_coloring(args.initial_report),
        colors=args.colors,
        solver_name=args.solver,
        max_rounds=args.max_rounds,
    )
    payload: dict[str, object] = {**stats, "coloring": coloring}
    _write_json(args.output, payload)
    print(json.dumps(stats, sort_keys=True))
    return 10 if coloring is not None else 0


if __name__ == "__main__":
    raise SystemExit(main())
