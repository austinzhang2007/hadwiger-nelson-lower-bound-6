"""CEGIS over exact five-color blocking points in the primary-D4 field."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
from typing import Callable, Sequence

import sympy as sp

from blocking_point_search import exact_five_color_blockers
from coloring_sat import ColoringSAT, parse_dimacs_edge, write_dimacs_edge
from d4_exact import (
    D4CoordinateBackend,
    D4ProjectivePoint,
    primary_d4_coordinate_backend,
)
from exact_geometry import (
    Point,
    parse_mathematica_vertices,
    parse_real_expr,
    serialize_graph,
)
from graph_verifier import load_graph_json


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _projective_payload(point: D4ProjectivePoint) -> dict[str, object]:
    def coefficients(value: Sequence[sp.Rational]) -> list[str]:
        return [sp.sstr(coefficient) for coefficient in value]

    return {
        "numerator_x_coefficients": coefficients(point.numerator_x),
        "numerator_y_coefficients": coefficients(point.numerator_y),
        "denominator_coefficients": coefficients(point.denominator),
        "basis_masks": list(range(32)),
    }


def _projective_from_payload(
    payload: dict[str, object],
) -> D4ProjectivePoint:
    def coefficients(key: str) -> tuple[sp.Rational, ...]:
        raw = payload[key]
        if not isinstance(raw, list) or len(raw) != 32:
            raise ValueError(f"invalid projective coefficient vector {key}")
        return tuple(sp.Rational(str(value)) for value in raw)

    return D4ProjectivePoint(
        numerator_x=coefficients("numerator_x_coefficients"),
        numerator_y=coefficients("numerator_y_coefficients"),
        denominator=coefficients("denominator_coefficients"),
        fifth_neighbor=int(payload["fifth_neighbor_zero_based"]),
    )


def _write_current_artifacts(
    *,
    output_directory: Path,
    base_points: Sequence[Point],
    base_edges: set[tuple[int, int]],
    candidate_points: Sequence[Point],
    current_edges: set[tuple[int, int]],
    iterations: Sequence[dict[str, object]],
    coloring: Sequence[int] | None,
    symmetry_clique: Sequence[int],
    status: str,
    complete_unit_distance_graph: bool,
) -> dict[str, object]:
    vertex_count = len(base_points) + len(candidate_points)
    stem = f"heule{vertex_count}-d4-blockers"
    graph_path = output_directory / f"{stem}.graph.json"
    edge_path = output_directory / f"{stem}.edge"
    cnf_path = output_directory / f"{stem}-5.cnf"
    report_path = output_directory / f"{stem}.json"
    points = [*base_points, *candidate_points]
    _write_json(
        graph_path,
        serialize_graph(
            "Complete primary-D4 graph plus exact five-color blockers",
            points,
            current_edges,
        ),
    )
    write_dimacs_edge(edge_path, vertex_count, current_edges)
    encoding = ColoringSAT(
        vertex_count,
        current_edges,
        5,
        break_color_symmetry=False,
        symmetry_clique=symmetry_clique,
    )
    encoding.write_dimacs_cnf(cnf_path)
    payload: dict[str, object] = {
        "status": status,
        "lower_bound_6_proved": False,
        "vertices": vertex_count,
        "edges": len(current_edges),
        "base_vertices": len(base_points),
        "base_edges": len(base_edges),
        "blocker_count": len(candidate_points),
        "complete_unit_distance_graph": complete_unit_distance_graph,
        "edge_completeness_scope": (
            "all base edges from the independently reconstructed stage3 graph; "
            + (
                "all blocker-to-base and blocker-to-blocker edges by exact "
                "32-dimensional primary-D4 algebra"
                if complete_unit_distance_graph
                else "five exact color-blocking generator edges per newly "
                "listed blocker; additional unit edges intentionally omitted"
            )
        ),
        "symmetry_clique_zero_based": list(symmetry_clique),
        "iterations": list(iterations),
        "whole_graph_5_sat": coloring is not None,
        "whole_graph_coloring": list(coloring) if coloring is not None else None,
        "graph_path": str(graph_path),
        "edge_path": str(edge_path),
        "cnf_path": str(cnf_path),
        "cnf_variables": encoding.variable_count,
        "cnf_clauses": len(encoding.clauses),
    }
    _write_json(report_path, payload)
    payload["sha256"] = {
        "graph": _sha256(graph_path),
        "edge": _sha256(edge_path),
        "cnf": _sha256(cnf_path),
    }
    _write_json(report_path, payload)
    return payload


def run_d4_blocker_cegis(
    *,
    graph_path: str | Path,
    coloring_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    output_directory: str | Path,
    max_iterations: int = 20,
    symmetry_clique: Sequence[int] = (548, 1149, 668),
    progress: Callable[[str, dict[str, object]], None] | None = None,
    resume_report_path: str | Path | None = None,
    blockers_per_round: int = 1,
    complete_edges_during_search: bool = True,
) -> dict[str, object]:
    """Run exact single-point-blocker CEGIS from a validated stage3 coloring."""

    started = time.monotonic()
    graph_path = Path(graph_path)
    output_directory = Path(output_directory)
    base_points, base_edges = load_graph_json(graph_path)
    with Path(coloring_path).open(encoding="utf-8") as handle:
        coloring_payload = json.load(handle)
    coloring = tuple(int(color) for color in coloring_payload["coloring"])
    if not ColoringSAT(
        len(base_points),
        base_edges,
        5,
        break_color_symmetry=False,
    ).validate(coloring):
        raise ValueError("input coloring is invalid")
    with Path(primary_report_path).open(encoding="utf-8") as handle:
        primary_report = json.load(handle)
    selected = [
        int(index)
        for index in primary_report["selected_source_r_indices_zero_based"]
    ]
    source_points = parse_mathematica_vertices(source_vertex_path)
    base_vertex_count = len(base_points) - len(selected) * len(source_points)
    backend = primary_d4_coordinate_backend(
        base_points[:base_vertex_count],
        source_points,
        selected,
    )
    if len(backend.coordinates) != len(base_points):
        raise AssertionError("D4 backend point count mismatch")

    current_edges = set(base_edges)
    candidate_points: list[Point] = []
    candidate_projective: list[D4ProjectivePoint] = []
    iterations: list[dict[str, object]] = []
    if resume_report_path is not None:
        with Path(resume_report_path).open(encoding="utf-8") as handle:
            resume_report = json.load(handle)
        if int(resume_report["base_vertices"]) != len(base_points):
            raise ValueError("resume report uses a different base graph")
        resume_vertex_count, resume_edges = parse_dimacs_edge(
            resume_report["edge_path"]
        )
        current_edges = set(resume_edges)
        iterations = list(resume_report["iterations"])
        candidate_projective = [
            backend.algebra.normalize_projective(
                _projective_from_payload(
                    {
                        **record["projective_certificate"],
                        "fifth_neighbor_zero_based": record[
                            "projective_certificate"
                        ].get(
                            "fifth_neighbor_zero_based",
                            record.get("generator_4_zero_based", -1),
                        ),
                    }
                )
            )
            for record in iterations
        ]
        candidate_points = [
            backend.materialize(projective)
            for projective in candidate_projective
        ]
        if resume_vertex_count != len(base_points) + len(candidate_points):
            raise ValueError("resume EDGE vertex count differs")
        if len(candidate_points) != len(candidate_projective):
            raise ValueError("resume candidate certificate count differs")
        coloring = tuple(
            int(color) for color in resume_report["whole_graph_coloring"]
        )
        if not ColoringSAT(
            resume_vertex_count,
            current_edges,
            5,
            break_color_symmetry=False,
        ).validate(coloring):
            raise ValueError("resume coloring is invalid")
    status = "ITERATION_LIMIT_WITH_5_COLORING"
    latest_payload: dict[str, object] = {}
    if blockers_per_round < 1:
        raise ValueError("blockers_per_round must be positive")
    projective_by_numeric_key: dict[tuple[float, float], list[int]] = {}
    for index, point in enumerate(candidate_points):
        approximate = point.approximate()
        key = (round(approximate[0], 10), round(approximate[1], 10))
        projective_by_numeric_key.setdefault(key, []).append(index)

    for round_index in range(1, max_iterations + 1):
        scan_started = time.monotonic()
        blockers, stats = exact_five_color_blockers(
            base_points,
            coloring[: len(base_points)],
            exact_backend=backend,
            max_blockers=blockers_per_round,
        )
        if not blockers:
            status = "SINGLE_POINT_POOL_STAGNATED_WITH_5_COLORING"
            break
        base_coloring = coloring[: len(base_points)]
        round_records: list[dict[str, object]] = []
        for blocker in blockers:
            neighbor_by_color = {
                int(base_coloring[index]): index for index in blocker.neighbors
            }
            second_generator = (
                neighbor_by_color[2],
                neighbor_by_color[3],
            )
            projective = backend.recover_projective(
                blocker.generator_pair,
                second_generator,
                (neighbor_by_color[4],),
            )
            if projective is None:
                raise AssertionError("saved blocker lost its exact D4 certificate")
            approximate = blocker.point.approximate()
            numeric_key = (
                round(approximate[0], 10),
                round(approximate[1], 10),
            )
            if any(
                backend.projective_equal(
                    candidate_projective[index],
                    projective,
                )
                for index in projective_by_numeric_key.get(numeric_key, ())
            ):
                continue
            complete_base_neighbors = (
                backend.unit_neighbors(projective)
                if complete_edges_during_search
                else tuple(sorted(blocker.neighbors))
            )
            if set(
                base_coloring[index] for index in complete_base_neighbors
            ) != set(range(5)):
                raise AssertionError("exact blocker no longer sees all five colors")
            if blocker.point in set(base_points):
                raise AssertionError("CEGIS selected a duplicate point")

            vertex = len(base_points) + len(candidate_points)
            current_edges.update(
                (min(index, vertex), max(index, vertex))
                for index in complete_base_neighbors
            )
            interaction_neighbors: list[int] = []
            for index, previous in enumerate(candidate_projective):
                if complete_edges_during_search and backend.projective_unit(
                    previous,
                    projective,
                ):
                    previous_vertex = len(base_points) + index
                    current_edges.add((previous_vertex, vertex))
                    interaction_neighbors.append(previous_vertex)
            candidate_points.append(blocker.point)
            candidate_projective.append(projective)
            projective_by_numeric_key.setdefault(numeric_key, []).append(
                len(candidate_projective) - 1
            )
            candidate_number = len(iterations) + len(round_records) + 1
            round_records.append(
                {
                    "iteration": candidate_number,
                    "cegis_round": round_index,
                    "candidate_vertex_zero_based": vertex,
                    "x_exact": sp.sstr(blocker.point.x),
                    "y_exact": sp.sstr(blocker.point.y),
                    "generator_01_zero_based": list(blocker.generator_pair),
                    "generator_23_zero_based": list(second_generator),
                    "generator_4_zero_based": neighbor_by_color[4],
                    "complete_base_neighbors_zero_based": list(
                        complete_base_neighbors
                    ),
                    "complete_base_degree": len(complete_base_neighbors),
                    "candidate_interaction_neighbors_zero_based": (
                        interaction_neighbors
                    ),
                    "projective_certificate": {
                        **_projective_payload(projective),
                        "fifth_neighbor_zero_based": (
                            projective.fifth_neighbor
                        ),
                    },
                }
            )
            if progress is not None:
                progress(
                    "blocker_added",
                    {
                        "iteration": candidate_number,
                        "round": round_index,
                        "vertex": vertex,
                        "base_degree": len(complete_base_neighbors),
                        "candidate_interactions": len(interaction_neighbors),
                    },
                )

        solve_started = time.monotonic()
        if not round_records:
            status = "LISTED_SINGLE_POINT_POOL_REPEATED_WITH_5_COLORING"
            break
        phases = [
            vertex_index * 5 + color + 1
            if color == coloring[vertex_index]
            else -(vertex_index * 5 + color + 1)
            for vertex_index in range(len(coloring))
            for color in range(5)
        ]
        next_coloring = ColoringSAT(
            len(base_points) + len(candidate_points),
            current_edges,
            5,
            break_color_symmetry=False,
            symmetry_clique=symmetry_clique,
        ).solve(phases=phases)
        solve_seconds = time.monotonic() - solve_started
        if progress is not None:
            progress(
                "sat_result",
                {
                    "round": round_index,
                    "candidate_count": len(candidate_points),
                    "five_sat": next_coloring is not None,
                    "solve_seconds": solve_seconds,
                },
            )
        for offset, record in enumerate(round_records):
            record["blocking_search"] = asdict(stats) if offset == 0 else None
            record["scan_seconds"] = (
                time.monotonic() - scan_started - solve_seconds
                if offset == 0
                else 0.0
            )
            record["solve_seconds"] = (
                solve_seconds if offset == len(round_records) - 1 else 0.0
            )
            record["five_sat_after_add"] = next_coloring is not None
        iterations.extend(round_records)
        if next_coloring is None:
            coloring = ()
            status = "EXACT_COMPLETE_5_UNSAT_WITHOUT_PROOF"
        else:
            coloring = next_coloring
            status = "RUNNING_WITH_5_COLORING"
        latest_payload = _write_current_artifacts(
            output_directory=output_directory,
            base_points=base_points,
            base_edges=base_edges,
            candidate_points=candidate_points,
            current_edges=current_edges,
            iterations=iterations,
            coloring=coloring if coloring else None,
            symmetry_clique=symmetry_clique,
            status=status,
            complete_unit_distance_graph=complete_edges_during_search,
        )
        if next_coloring is None:
            break

    if status == "RUNNING_WITH_5_COLORING":
        status = "ITERATION_LIMIT_WITH_5_COLORING"
        latest_payload = _write_current_artifacts(
            output_directory=output_directory,
            base_points=base_points,
            base_edges=base_edges,
            candidate_points=candidate_points,
            current_edges=current_edges,
            iterations=iterations,
            coloring=coloring,
            symmetry_clique=symmetry_clique,
            status=status,
            complete_unit_distance_graph=complete_edges_during_search,
        )
    if status in {
        "SINGLE_POINT_POOL_STAGNATED_WITH_5_COLORING",
        "LISTED_SINGLE_POINT_POOL_REPEATED_WITH_5_COLORING",
        "ITERATION_LIMIT_WITH_5_COLORING",
    }:
        latest_payload = _write_current_artifacts(
            output_directory=output_directory,
            base_points=base_points,
            base_edges=base_edges,
            candidate_points=candidate_points,
            current_edges=current_edges,
            iterations=iterations,
            coloring=coloring,
            symmetry_clique=symmetry_clique,
            status=status,
            complete_unit_distance_graph=complete_edges_during_search,
        )
    latest_payload["elapsed_seconds"] = time.monotonic() - started
    report_path = Path(str(latest_payload["graph_path"])).with_suffix("").with_suffix(
        ".json"
    )
    _write_json(report_path, latest_payload)
    return latest_payload
