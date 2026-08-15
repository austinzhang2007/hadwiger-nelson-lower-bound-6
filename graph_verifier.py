"""Independent finite unit-distance graph and DRAT certificate verifier."""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from exact_geometry import Point, complete_unit_edges, parse_mathematica_vertices
from exact_geometry import parse_real_expr


Edge = tuple[int, int]


@dataclass(frozen=True)
class GeometryReport:
    vertex_count: int
    edge_count: int
    duplicate_vertices: tuple[tuple[int, int], ...]
    missing_edges: tuple[Edge, ...]
    extra_edges: tuple[Edge, ...]

    @property
    def valid(self) -> bool:
        return not (
            self.duplicate_vertices or self.missing_edges or self.extra_edges
        )


@dataclass(frozen=True)
class ProofReport:
    valid: bool
    returncode: int
    summary: str


def _parse_dimacs_edge(path: str | Path) -> tuple[int, set[Edge]]:
    vertex_count = declared_edge_count = None
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
                vertex_count = int(fields[2])
                declared_edge_count = int(fields[3])
            elif fields[0] == "e":
                if vertex_count is None or len(fields) != 3:
                    raise ValueError(f"{path}:{line_number}: malformed edge")
                u, v = int(fields[1]) - 1, int(fields[2]) - 1
                if u == v or not (0 <= u < vertex_count and 0 <= v < vertex_count):
                    raise ValueError(f"{path}:{line_number}: invalid endpoint")
                edge = (min(u, v), max(u, v))
                if edge in edges:
                    raise ValueError(f"{path}:{line_number}: duplicate edge")
                edges.add(edge)
            else:
                raise ValueError(f"{path}:{line_number}: unknown record")
    if vertex_count is None or declared_edge_count is None:
        raise ValueError(f"{path}: missing p edge header")
    if len(edges) != declared_edge_count:
        raise ValueError(
            f"{path}: declares {declared_edge_count} edges, has {len(edges)}"
        )
    return vertex_count, edges


def _duplicate_vertices(points: list[Point]) -> tuple[tuple[int, int], ...]:
    first_occurrence: dict[Point, int] = {}
    duplicates: list[tuple[int, int]] = []
    for index, point in enumerate(points):
        if point in first_occurrence:
            duplicates.append((first_occurrence[point], index))
        else:
            first_occurrence[point] = index
    return tuple(duplicates)


def verify_points_and_edges(
    points: list[Point], declared_edges: Iterable[Edge]
) -> GeometryReport:
    normalized = {(min(u, v), max(u, v)) for u, v in declared_edges}
    if any(
        u == v or u < 0 or v >= len(points)
        for u, v in normalized
    ):
        raise ValueError("declared edge endpoint outside graph")
    exact_edges = complete_unit_edges(points)
    return GeometryReport(
        vertex_count=len(points),
        edge_count=len(exact_edges),
        duplicate_vertices=_duplicate_vertices(points),
        missing_edges=tuple(sorted(exact_edges - normalized)),
        extra_edges=tuple(sorted(normalized - exact_edges)),
    )


def verify_author_graph(
    coordinate_path: str | Path, edge_path: str | Path
) -> GeometryReport:
    points = parse_mathematica_vertices(coordinate_path)
    vertex_count, edges = _parse_dimacs_edge(edge_path)
    if len(points) != vertex_count:
        raise ValueError(
            f"coordinate count {len(points)} != edge header {vertex_count}"
        )
    return verify_points_and_edges(points, edges)


def load_graph_json(path: str | Path) -> tuple[list[Point], set[Edge]]:
    with Path(path).open(encoding="utf-8") as handle:
        data = json.load(handle)
    if data.get("schema") != "udg-exact-v1":
        raise ValueError("unsupported graph schema")
    if data.get("coordinate_format") != "sympy-radical-v1":
        raise ValueError("unsupported coordinate format")
    vertices = data.get("vertices")
    if not isinstance(vertices, list):
        raise ValueError("vertices must be a list")
    points: list[Point] = []
    for expected_id, vertex in enumerate(vertices, 1):
        if not isinstance(vertex, dict) or vertex.get("id") != expected_id:
            raise ValueError("vertex ids must be consecutive and one-based")
        points.append(
            Point(
                parse_real_expr(str(vertex["x"])),
                parse_real_expr(str(vertex["y"])),
            )
        )
    raw_edges = data.get("edges")
    if not isinstance(raw_edges, list):
        raise ValueError("edges must be a list")
    edges: set[Edge] = set()
    for raw_edge in raw_edges:
        if (
            not isinstance(raw_edge, list)
            or len(raw_edge) != 2
            or not all(isinstance(value, int) for value in raw_edge)
        ):
            raise ValueError("each edge must be a two-integer list")
        u, v = raw_edge[0] - 1, raw_edge[1] - 1
        edge = (min(u, v), max(u, v))
        if edge in edges:
            raise ValueError("duplicate edge")
        edges.add(edge)
    return points, edges


def verify_graph_json(path: str | Path) -> GeometryReport:
    points, edges = load_graph_json(path)
    return verify_points_and_edges(points, edges)


def verify_drat(
    checker_path: str | Path,
    cnf_path: str | Path,
    proof_path: str | Path,
) -> ProofReport:
    process = subprocess.run(
        [str(Path(checker_path)), str(Path(cnf_path)), str(Path(proof_path))],
        check=False,
        capture_output=True,
        text=True,
    )
    output = process.stdout + process.stderr
    valid = process.returncode == 0 and "VERIFIED" in output
    summary_lines = [
        line.strip()
        for line in output.splitlines()
        if "VERIFIED" in line or "NOT VERIFIED" in line or "s VERIFIED" in line
    ]
    return ProofReport(valid, process.returncode, " | ".join(summary_lines))


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--graph-json", type=Path)
    source.add_argument("--coordinates", type=Path)
    parser.add_argument("--edge-file", type=Path)
    parser.add_argument("--drat-trim", type=Path)
    parser.add_argument("--cnf", type=Path)
    parser.add_argument("--drat", type=Path)
    args = parser.parse_args()

    if args.coordinates and not args.edge_file:
        parser.error("--coordinates requires --edge-file")
    geometry = (
        verify_graph_json(args.graph_json)
        if args.graph_json
        else verify_author_graph(args.coordinates, args.edge_file)
    )
    result: dict[str, object] = {"geometry": asdict(geometry)}
    valid = geometry.valid
    proof_arguments = (args.drat_trim, args.cnf, args.drat)
    if any(proof_arguments):
        if not all(proof_arguments):
            parser.error("--drat-trim, --cnf and --drat must be supplied together")
        proof = verify_drat(*proof_arguments)
        result["proof"] = asdict(proof)
        valid = valid and proof.valid
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if valid else 1


if __name__ == "__main__":
    raise SystemExit(_main())

