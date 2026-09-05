"""Lean equisatisfiable graph-coloring CNF without at-most-one clauses.

The clauses require every vertex to have at least one color and adjacent
vertices to share no color.  A satisfying model may assign several colors to
one vertex, but choosing any one true color per vertex always yields an
ordinary proper coloring.  Thus the encoding is equisatisfiable with the
standard exactly-one encoding and is stronger evidence when it is UNSAT.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

from coloring_sat import Edge, parse_dimacs_edge


def _variable(vertex: int, color: int, colors: int) -> int:
    return vertex * colors + color + 1


def _normalized_graph(
    vertex_count: int,
    edges: Iterable[Edge],
    colors: int,
    symmetry_clique: Sequence[int],
) -> set[Edge]:
    if vertex_count < 0 or colors < 1:
        raise ValueError("invalid vertex or color count")
    normalized = {(min(left, right), max(left, right)) for left, right in edges}
    if any(
        left < 0 or right >= vertex_count or left == right
        for left, right in normalized
    ):
        raise ValueError("edge endpoint outside graph")
    if (
        len(symmetry_clique) > colors
        or len(set(symmetry_clique)) != len(symmetry_clique)
        or any(vertex < 0 or vertex >= vertex_count for vertex in symmetry_clique)
    ):
        raise ValueError("invalid symmetry clique")
    if any(
        (min(left, right), max(left, right)) not in normalized
        for left, right in combinations(symmetry_clique, 2)
    ):
        raise ValueError("symmetry_clique must induce a clique")
    return normalized


def lean_coloring_clauses(
    vertex_count: int,
    edges: Iterable[Edge],
    colors: int,
    *,
    symmetry_clique: Sequence[int] = (),
) -> list[list[int]]:
    """Build the lean coloring clauses in deterministic order."""

    normalized = _normalized_graph(
        vertex_count, edges, colors, symmetry_clique
    )
    clauses = [
        [_variable(vertex, color, colors) for color in range(colors)]
        for vertex in range(vertex_count)
    ]
    clauses.extend(
        [
            -_variable(left, color, colors),
            -_variable(right, color, colors),
        ]
        for left, right in sorted(normalized)
        for color in range(colors)
    )
    clauses.extend(
        [_variable(vertex, color, colors)]
        for color, vertex in enumerate(symmetry_clique)
    )
    return clauses


def decode_lean_coloring_model(
    model: Sequence[int],
    *,
    vertex_count: int,
    edges: Iterable[Edge],
    colors: int,
) -> tuple[int, ...]:
    """Choose one true color per vertex and validate every graph edge."""

    positive = {literal for literal in model if literal > 0}
    coloring: list[int] = []
    for vertex in range(vertex_count):
        selected = [
            color
            for color in range(colors)
            if _variable(vertex, color, colors) in positive
        ]
        if not selected:
            raise ValueError(f"vertex {vertex} has no selected color")
        coloring.append(selected[0])
    normalized = {(min(left, right), max(left, right)) for left, right in edges}
    if any(coloring[left] == coloring[right] for left, right in normalized):
        raise ValueError("decoded lean model contains a monochromatic edge")
    return tuple(coloring)


def write_lean_coloring_cnf(
    edge_path: str | Path,
    target_path: str | Path,
    *,
    colors: int = 5,
    symmetry_clique: Sequence[int] = (),
) -> dict[str, int]:
    """Stream a lean coloring DIMACS file from a DIMACS edge graph."""

    vertex_count, edges = parse_dimacs_edge(edge_path)
    normalized = _normalized_graph(
        vertex_count, edges, colors, symmetry_clique
    )
    variables = vertex_count * colors
    clause_count = (
        vertex_count + len(normalized) * colors + len(symmetry_clique)
    )
    target = Path(target_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    with temporary.open("w", encoding="ascii") as handle:
        handle.write(f"p cnf {variables} {clause_count}\n")
        for vertex in range(vertex_count):
            handle.write(
                " ".join(
                    str(_variable(vertex, color, colors))
                    for color in range(colors)
                )
            )
            handle.write(" 0\n")
        for left, right in sorted(normalized):
            for color in range(colors):
                handle.write(
                    f"-{_variable(left, color, colors)} "
                    f"-{_variable(right, color, colors)} 0\n"
                )
        for color, vertex in enumerate(symmetry_clique):
            handle.write(f"{_variable(vertex, color, colors)} 0\n")
    temporary.replace(target)
    return {
        "vertices": vertex_count,
        "edges": len(normalized),
        "colors": colors,
        "variables": variables,
        "clauses": clause_count,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--colors", type=int, default=5)
    parser.add_argument("--symmetry-clique", type=int, nargs="*", default=[])
    args = parser.parse_args()
    stats = write_lean_coloring_cnf(
        args.edge,
        args.output,
        colors=args.colors,
        symmetry_clique=args.symmetry_clique,
    )
    print(stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
