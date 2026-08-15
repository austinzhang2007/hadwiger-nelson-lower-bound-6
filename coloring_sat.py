"""Incremental SAT encoding for exact graph coloring.

Vertex and color indices are zero-based in Python.  DIMACS variables are
``vertex * color_count + color + 1``.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from pysat.solvers import Solver


Edge = tuple[int, int]
Coloring = tuple[int, ...]


def parse_dimacs_edge(path: str | Path) -> tuple[int, set[Edge]]:
    vertex_count = edge_count = None
    edges: set[Edge] = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("c"):
                continue
            fields = line.split()
            if fields[0] == "p":
                if len(fields) != 4 or fields[1] != "edge":
                    raise ValueError(f"{path}:{line_number}: malformed header")
                vertex_count, edge_count = int(fields[2]), int(fields[3])
            elif fields[0] == "e":
                if vertex_count is None or len(fields) != 3:
                    raise ValueError(f"{path}:{line_number}: malformed edge")
                u, v = int(fields[1]) - 1, int(fields[2]) - 1
                if not (0 <= u < vertex_count and 0 <= v < vertex_count) or u == v:
                    raise ValueError(f"{path}:{line_number}: invalid edge")
                edges.add((min(u, v), max(u, v)))
            else:
                raise ValueError(f"{path}:{line_number}: unknown record")
    if vertex_count is None or edge_count is None:
        raise ValueError(f"{path}: missing p edge header")
    if len(edges) != edge_count:
        raise ValueError(
            f"{path}: header declares {edge_count} edges, parsed {len(edges)}"
        )
    return vertex_count, edges


def write_dimacs_edge(
    path: str | Path, vertex_count: int, edges: Iterable[Edge]
) -> None:
    normalized = sorted({(min(u, v), max(u, v)) for u, v in edges})
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write(f"p edge {vertex_count} {len(normalized)}\n")
        for u, v in normalized:
            handle.write(f"e {u + 1} {v + 1}\n")


class ColoringSAT:
    """Exactly-one graph coloring with full color-permutation symmetry breaking."""

    def __init__(
        self,
        vertex_count: int,
        edges: Iterable[Edge],
        color_count: int,
        *,
        break_color_symmetry: bool = True,
        symmetry_clique: Sequence[int] = (),
        solver_name: str = "cadical195",
    ) -> None:
        if vertex_count < 0 or color_count <= 0:
            raise ValueError("invalid vertex or color count")
        self.vertex_count = vertex_count
        self.color_count = color_count
        raw_edges = list(edges)
        if any(
            u < 0
            or v < 0
            or u >= vertex_count
            or v >= vertex_count
            or u == v
            for u, v in raw_edges
        ):
            raise ValueError("edge endpoint outside graph")
        self.edges = {
            (min(u, v), max(u, v)) for u, v in raw_edges
        }
        self.break_color_symmetry = break_color_symmetry
        self.symmetry_clique = tuple(symmetry_clique)
        if self.break_color_symmetry and self.symmetry_clique:
            raise ValueError(
                "choose either color precedence or a symmetry clique"
            )
        if (
            len(self.symmetry_clique) > self.color_count
            or len(set(self.symmetry_clique)) != len(self.symmetry_clique)
            or any(
                vertex < 0 or vertex >= self.vertex_count
                for vertex in self.symmetry_clique
            )
        ):
            raise ValueError("invalid symmetry clique")
        if any(
            (min(left, right), max(left, right)) not in self.edges
            for left, right in combinations(self.symmetry_clique, 2)
        ):
            raise ValueError("symmetry_clique must induce a clique")
        self.solver_name = solver_name
        self.clauses = self._build_clauses()

    def variable(self, vertex: int, color: int) -> int:
        return vertex * self.color_count + color + 1

    @property
    def variable_count(self) -> int:
        return self.vertex_count * self.color_count

    def _build_clauses(self) -> list[list[int]]:
        clauses: list[list[int]] = []
        for vertex in range(self.vertex_count):
            variables = [
                self.variable(vertex, color)
                for color in range(self.color_count)
            ]
            clauses.append(variables)
            clauses.extend([-left, -right] for left, right in combinations(variables, 2))
        for u, v in sorted(self.edges):
            for color in range(self.color_count):
                clauses.append(
                    [-self.variable(u, color), -self.variable(v, color)]
                )
        if self.break_color_symmetry and self.vertex_count:
            # Restricted-growth-string precedence: the first use of color c
            # occurs strictly after the first use of color c-1.
            for color in range(1, self.color_count):
                preceding: list[int] = []
                for vertex in range(self.vertex_count):
                    clauses.append(
                        [-self.variable(vertex, color), *preceding]
                    )
                    preceding.append(self.variable(vertex, color - 1))
        for color, vertex in enumerate(self.symmetry_clique):
            clauses.append([self.variable(vertex, color)])
        return clauses

    def _decode(self, model: Sequence[int]) -> Coloring:
        positive = {literal for literal in model if literal > 0}
        coloring: list[int] = []
        for vertex in range(self.vertex_count):
            selected = [
                color
                for color in range(self.color_count)
                if self.variable(vertex, color) in positive
            ]
            if len(selected) != 1:
                raise RuntimeError(
                    f"SAT model selects {len(selected)} colors for vertex {vertex}"
                )
            coloring.append(selected[0])
        result = tuple(coloring)
        if not self.validate(result):
            raise RuntimeError("solver returned an invalid coloring")
        return result

    def validate(self, coloring: Sequence[int]) -> bool:
        return (
            len(coloring) == self.vertex_count
            and all(0 <= color < self.color_count for color in coloring)
            and all(coloring[u] != coloring[v] for u, v in self.edges)
        )

    def solve(
        self,
        assumptions: Sequence[int] = (),
        *,
        phases: Sequence[int] = (),
    ) -> Coloring | None:
        with Solver(
            name=self.solver_name, bootstrap_with=self.clauses
        ) as solver:
            if phases:
                solver.set_phases(list(phases))
            if not solver.solve(assumptions=list(assumptions)):
                return None
            return self._decode(solver.get_model())

    def enumerate_colorings(self, limit: int | None = None) -> Iterator[Coloring]:
        if limit is not None and limit < 0:
            raise ValueError("limit must be nonnegative")
        count = 0
        with Solver(
            name=self.solver_name, bootstrap_with=self.clauses
        ) as solver:
            while (limit is None or count < limit) and solver.solve():
                coloring = self._decode(solver.get_model())
                yield coloring
                count += 1
                solver.add_clause(
                    [
                        -self.variable(vertex, color)
                        for vertex, color in enumerate(coloring)
                    ]
                )

    def write_dimacs_cnf(self, path: str | Path) -> None:
        with Path(path).open("w", encoding="ascii") as handle:
            handle.write(
                f"p cnf {self.variable_count} {len(self.clauses)}\n"
            )
            for clause in self.clauses:
                handle.write(" ".join(map(str, clause)))
                handle.write(" 0\n")

    def solve_with_drat(
        self,
        proof_path: str | Path,
        *,
        proof_solver: str = "glucose4",
    ) -> Coloring | None:
        """Solve and, only for UNSAT, write the solver's textual DRUP/DRAT trace."""

        with Solver(
            name=proof_solver,
            bootstrap_with=self.clauses,
            with_proof=True,
        ) as solver:
            if solver.solve():
                return self._decode(solver.get_model())
            proof = solver.get_proof()
        if proof is None:
            raise RuntimeError(f"{proof_solver} did not return a proof trace")
        with Path(proof_path).open("w", encoding="ascii") as handle:
            for line in proof:
                handle.write(line.rstrip())
                handle.write("\n")
        return None


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("edge_file", type=Path)
    parser.add_argument("-k", "--colors", type=int, required=True)
    parser.add_argument("--enumerate", type=int, default=0)
    parser.add_argument("--cnf", type=Path)
    parser.add_argument("--proof", type=Path)
    args = parser.parse_args()

    vertex_count, edges = parse_dimacs_edge(args.edge_file)
    encoding = ColoringSAT(vertex_count, edges, args.colors)
    if args.cnf:
        encoding.write_dimacs_cnf(args.cnf)
    if args.enumerate:
        models = list(encoding.enumerate_colorings(args.enumerate))
        print(f"models={len(models)}")
        return 0
    coloring = (
        encoding.solve_with_drat(args.proof) if args.proof else encoding.solve()
    )
    print("UNSAT" if coloring is None else "SAT")
    return 20 if coloring is None else 10


if __name__ == "__main__":
    raise SystemExit(_main())
