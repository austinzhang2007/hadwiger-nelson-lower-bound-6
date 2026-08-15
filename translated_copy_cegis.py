"""CEGIS with exact quadratic translations of the Heule 529-point graph.

The search-facing construction deliberately avoids a generic scan in the
degree-16 compositum.  For the certified non-base quadratic translation used
here, every base-to-copy unit edge has one of two exact anchors; see
``docs/plans/2026-07-30-translated-copy-cegis-design.md``.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp

from cegis import candidate_extension
from coloring_sat import ColoringSAT, parse_dimacs_edge
from exact_geometry import Point, parse_mathematica_vertices, serialize_graph
from graph_verifier import load_graph_json


Edge = tuple[int, int]


@dataclass(frozen=True, slots=True)
class CircleTranslation:
    """One intersection of the unit circles centered at 0 and ``r``."""

    point: Point
    h_squared: sp.Expr
    minimal_polynomial: sp.Poly
    galois_group_is_nonabelian: bool
    branch: int


@dataclass(frozen=True, slots=True)
class TranslationCandidate:
    """One exact base-list constraint induced by a source translation."""

    source_r_index: int
    r: Point
    second_neighbor_count: int
    neighbors: dict[int, tuple[int, ...]]


def circle_translation(r: Point, *, branch: int) -> CircleTranslation:
    """Construct and certify the chosen quadratic circle translation.

    A nonabelian Galois group for the minimal polynomial of ``h`` proves
    ``h`` is outside every multiquadratic (hence abelian Galois) base field.
    """

    if branch not in (0, 1):
        raise ValueError("branch must be 0 or 1")
    d2 = sp.radsimp(sp.expand(r.x * r.x + r.y * r.y))
    if d2 == 0 or sp.N(d2) >= 4:
        raise ValueError("circle centers must have exact distance in (0, 2)")
    h_squared = sp.radsimp(sp.cancel((4 - d2) / (4 * d2)))
    h = sp.sqrt(h_squared)
    sign = 1 if branch == 1 else -1
    point = Point(
        r.x / 2 + sign * r.y * h,
        r.y / 2 - sign * r.x * h,
    )
    variable = sp.symbols("x")
    polynomial = sp.Poly(
        sp.minimal_polynomial(h, variable),
        variable,
        domain=sp.QQ,
    )
    if polynomial.degree() > 6:
        raise ValueError(
            "SymPy's exact Galois-group certificate supports degree at most 6"
        )
    group, _ = polynomial.galois_group()
    nonabelian = not bool(group.is_abelian)
    if not nonabelian:
        raise ValueError(
            "translation has no nonabelian certificate outside the base field"
        )
    return CircleTranslation(
        point=point,
        h_squared=h_squared,
        minimal_polynomial=polynomial,
        galois_group_is_nonabelian=nonabelian,
        branch=branch,
    )


@lru_cache(maxsize=None)
def _is_in_q_sqrt_3_5_11(expression: sp.Expr) -> bool:
    """Recognize the explicit multiquadratic coordinate representation."""

    s3, s5, s11 = sp.symbols("s3 s5 s11")
    replacements = {
        sp.sqrt(165): s3 * s5 * s11,
        sp.sqrt(55): s5 * s11,
        sp.sqrt(33): s3 * s11,
        sp.sqrt(15): s3 * s5,
        sp.sqrt(11): s11,
        sp.sqrt(5): s5,
        sp.sqrt(3): s3,
    }
    replaced = sp.expand(expression.xreplace(replacements))
    if replaced.has(sp.sqrt):
        return False
    try:
        polynomial = sp.Poly(replaced, s3, s5, s11, domain=sp.QQ)
    except (sp.PolynomialError, sp.CoercionFailed):
        return False
    return all(all(exponent <= 1 for exponent in term) for term in polynomial.monoms())


def _points_are_in_base_field(points: Sequence[Point]) -> bool:
    return all(
        _is_in_q_sqrt_3_5_11(coordinate)
        for point in points
        for coordinate in (point.x, point.y)
    )


def translated_points(
    source_points: Sequence[Point],
    translation: Point,
) -> list[Point]:
    """Translate exact source coordinates without numeric approximation."""

    return [
        Point(point.x + translation.x, point.y + translation.y)
        for point in source_points
    ]


def primary_d4_candidate_pool(
    base_points: Sequence[Point],
    source_points: Sequence[Point],
) -> list[TranslationCandidate]:
    """Return the 24 primary translations with the certified D4 extension.

    Their two squared lengths are algebraic conjugates.  Both produce
    ``8*x**4 - 50*x**2 + 41`` for the auxiliary square root.
    """

    lookup = {point: index for index, point in enumerate(base_points)}
    if len(lookup) != len(base_points):
        raise ValueError("base graph contains duplicate exact points")
    targets = (
        sp.Rational(1, 2) - sp.sqrt(33) / 18,
        sp.Rational(1, 2) + sp.sqrt(33) / 18,
    )
    candidates: list[TranslationCandidate] = []
    for r_index, r in enumerate(source_points):
        d2 = sp.radsimp(sp.expand(r.x * r.x + r.y * r.y))
        if not any(sp.radsimp(d2 - target) == 0 for target in targets):
            continue
        neighbors: dict[int, tuple[int, ...]] = {}
        second_count = 0
        for source_index, point in enumerate(source_points):
            first = lookup.get(point)
            if first is None:
                raise ValueError("the source graph is not contained in the base")
            indices = [first]
            second = lookup.get(Point(point.x + r.x, point.y + r.y))
            if second is not None and second != first:
                indices.append(second)
                second_count += 1
            neighbors[source_index] = tuple(indices)
        candidates.append(
            TranslationCandidate(
                source_r_index=r_index,
                r=r,
                second_neighbor_count=second_count,
                neighbors=neighbors,
            )
        )
    return sorted(
        candidates,
        key=lambda candidate: (
            candidate.second_neighbor_count,
            candidate.source_r_index,
        ),
        reverse=True,
    )


def enumerate_structured_cross_edges(
    base_points: Sequence[Point],
    source_points: Sequence[Point],
    r: Point,
    shifted_points: Sequence[Point],
    translation: CircleTranslation,
) -> set[Edge]:
    """Enumerate every base-to-copy unit edge by the two-anchor theorem."""

    if len(source_points) != len(shifted_points):
        raise ValueError("source and shifted point counts differ")
    if not translation.galois_group_is_nonabelian:
        raise ValueError("translation lacks an outside-base-field certificate")
    if not _points_are_in_base_field([*base_points, *source_points, r]):
        raise ValueError("the two-anchor theorem requires coordinates in K")
    expected_shifted = translated_points(source_points, translation.point)
    if list(shifted_points) != expected_shifted:
        raise ValueError("shifted coordinates do not use the certified translation")

    lookup: dict[Point, int] = {}
    for index, point in enumerate(base_points):
        if point in lookup:
            raise ValueError("base graph contains duplicate exact points")
        lookup[point] = index

    cross_edges: set[Edge] = set()
    for source_index, point in enumerate(source_points):
        for anchor in (point, Point(point.x + r.x, point.y + r.y)):
            base_index = lookup.get(anchor)
            if base_index is not None:
                cross_edges.add((base_index, source_index))

    for base_index, source_index in cross_edges:
        if (
            base_points[base_index].squared_distance(
                shifted_points[source_index]
            )
            != 1
        ):
            raise AssertionError("the structured enumeration produced a non-unit edge")
    return cross_edges


def certified_unit_shift_copy_edges(
    source_points: Sequence[Point],
    selected_r_indices: Sequence[int],
) -> tuple[set[tuple[int, int, int, int]], dict[tuple[int, int], int]]:
    """Certify interactions when two copy translations differ by one unit.

    For equal ``h_squared``, write ``s = r_b-r_a``.  The translation difference
    is ``s/2 + (s_y, -s_x)h``.  If its squared length is one, separating the
    coefficients of ``1,h`` shows that a cross-copy edge exists exactly when
    ``p_b=p_a`` or ``p_b=p_a-s``.
    """

    source_lookup = {point: index for index, point in enumerate(source_points)}
    if len(source_lookup) != len(source_points):
        raise ValueError("source graph contains duplicate exact points")
    translations = [
        circle_translation(source_points[index], branch=1)
        for index in selected_r_indices
    ]
    interactions: set[tuple[int, int, int, int]] = set()
    pair_counts: dict[tuple[int, int], int] = {}
    for left_copy in range(len(selected_r_indices)):
        for right_copy in range(left_copy + 1, len(selected_r_indices)):
            left_translation = translations[left_copy]
            right_translation = translations[right_copy]
            if (
                sp.radsimp(
                    left_translation.h_squared - right_translation.h_squared
                )
                != 0
            ):
                continue
            delta = Point(
                right_translation.point.x - left_translation.point.x,
                right_translation.point.y - left_translation.point.y,
            )
            if delta.squared_distance(Point(0, 0)) != 1:
                continue
            left_r = source_points[selected_r_indices[left_copy]]
            right_r = source_points[selected_r_indices[right_copy]]
            s = Point(right_r.x - left_r.x, right_r.y - left_r.y)
            pair_edges: set[tuple[int, int, int, int]] = set()
            for left_index, point in enumerate(source_points):
                right_indices = {
                    left_index,
                    source_lookup.get(Point(point.x - s.x, point.y - s.y)),
                }
                for right_index in right_indices:
                    if right_index is None:
                        continue
                    left_point = Point(
                        point.x + left_translation.point.x,
                        point.y + left_translation.point.y,
                    )
                    right_source = source_points[right_index]
                    right_point = Point(
                        right_source.x + right_translation.point.x,
                        right_source.y + right_translation.point.y,
                    )
                    if left_point.squared_distance(right_point) != 1:
                        raise AssertionError(
                            "unit-shift theorem produced a non-unit interaction"
                        )
                    pair_edges.add(
                        (left_copy, left_index, right_copy, right_index)
                    )
            interactions.update(pair_edges)
            pair_counts[(left_copy, right_copy)] = len(pair_edges)
    return interactions, pair_counts


def augmented_edges(
    base_vertex_count: int,
    base_edges: Iterable[Edge],
    source_edges: Iterable[Edge],
    cross_edges: Iterable[Edge],
) -> set[Edge]:
    """Build the listed-edge graph containing one translated source copy."""

    offset = base_vertex_count
    edges = set(base_edges)
    edges.update((offset + left, offset + right) for left, right in source_edges)
    edges.update((base, offset + source) for base, source in cross_edges)
    return edges


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def run_single_copy(
    *,
    base_graph_path: Path,
    base_report_path: Path,
    source_vertex_path: Path,
    source_edge_path: Path,
    output_graph_path: Path,
    output_report_path: Path,
    r_vertex: int = 98,
    branch: int = 1,
    color_count: int = 5,
) -> dict[str, object]:
    """Run one exact coloring→blocker→augmentation→new-coloring cycle."""

    started = time.monotonic()
    base_points, base_edges = load_graph_json(base_graph_path)
    with base_report_path.open(encoding="utf-8") as handle:
        base_report = json.load(handle)
    old_coloring = tuple(
        int(color)
        for color in base_report["whole_pool_failure_certificate"]["coloring"]
    )
    source_points = parse_mathematica_vertices(source_vertex_path)
    source_count, source_edges = parse_dimacs_edge(source_edge_path)
    if source_count != len(source_points):
        raise ValueError("source coordinate and edge counts disagree")

    r = base_points[r_vertex - 1]
    translation = circle_translation(r, branch=branch)
    shifted = translated_points(source_points, translation.point)
    if set(shifted) & set(base_points):
        raise ValueError("translated copy overlaps the base graph")
    cross = enumerate_structured_cross_edges(
        base_points,
        source_points,
        r,
        shifted,
        translation,
    )
    neighbors = {index: [] for index in range(source_count)}
    for base_index, source_index in cross:
        neighbors[source_index].append(base_index)
    fixed_extension = candidate_extension(
        neighbors,
        source_edges,
        old_coloring,
        color_count,
    )

    full_points = [*base_points, *shifted]
    full_edges = augmented_edges(
        len(base_points),
        base_edges,
        source_edges,
        cross,
    )
    solver = ColoringSAT(len(full_points), full_edges, color_count)
    coloring = solver.solve()
    if coloring is None:
        status = "LISTED_EDGE_GRAPH_5_UNSAT_WITHOUT_PROOF"
    else:
        status = "LISTED_EDGE_GRAPH_5_SAT"
        if not solver.validate(coloring):
            raise AssertionError("SAT solver returned an invalid coloring")

    _write_json(
        output_graph_path,
        serialize_graph(
            "2510-point closure plus one exact translated Heule-529 copy",
            full_points,
            full_edges,
        ),
    )
    report: dict[str, object] = {
        "status": status,
        "lower_bound_6_proved": False,
        "base_graph_path": str(base_graph_path),
        "source_vertex_path": str(source_vertex_path),
        "source_edge_path": str(source_edge_path),
        "r_vertex_one_based": r_vertex,
        "r": {"x": sp.sstr(r.x), "y": sp.sstr(r.y)},
        "translation_branch": branch,
        "translation": {
            "x": sp.sstr(translation.point.x),
            "y": sp.sstr(translation.point.y),
            "h_squared": sp.sstr(translation.h_squared),
            "minimal_polynomial": sp.sstr(
                translation.minimal_polynomial.as_expr()
            ),
            "galois_group_nonabelian": True,
        },
        "vertices": len(full_points),
        "listed_exact_unit_edges": len(full_edges),
        "cross_edges": len(cross),
        "two_base_neighbors": sum(
            len(indices) == 2 for indices in neighbors.values()
        ),
        "blocks_saved_base_coloring": fixed_extension is None,
        "whole_graph_5_sat": coloring is not None,
        "whole_graph_coloring": list(coloring) if coloring is not None else None,
        "output_graph_path": str(output_graph_path),
        "elapsed_seconds": time.monotonic() - started,
        "edge_completeness_scope": (
            "complete within base, within copy, and base-to-copy by the "
            "certified two-anchor theorem"
        ),
    }
    _write_json(output_report_path, report)
    return report


def run_primary_d4_cegis(
    *,
    base_graph_path: Path,
    base_report_path: Path,
    source_vertex_path: Path,
    source_edge_path: Path,
    output_graph_path: Path,
    output_report_path: Path,
    max_iterations: int = 24,
    color_count: int = 5,
) -> dict[str, object]:
    """Run CEGIS over the reproducible 24-member primary D4 pool."""

    started = time.monotonic()
    base_points, base_edges = load_graph_json(base_graph_path)
    with base_report_path.open(encoding="utf-8") as handle:
        base_report = json.load(handle)
    coloring = tuple(
        int(color)
        for color in base_report["whole_pool_failure_certificate"]["coloring"]
    )
    source_points = parse_mathematica_vertices(source_vertex_path)
    source_count, source_edges = parse_dimacs_edge(source_edge_path)
    if source_count != len(source_points):
        raise ValueError("source coordinate and edge counts disagree")

    pool = primary_d4_candidate_pool(base_points, source_points)
    points = list(base_points)
    point_indices = {point: index for index, point in enumerate(points)}
    edges = set(base_edges)
    selected: list[int] = []
    iterations: list[dict[str, object]] = []
    exhausted = False
    unsat = False

    for iteration_index in range(1, max_iterations + 1):
        scan_started = time.monotonic()
        blockers = [
            candidate
            for candidate in pool
            if candidate.source_r_index not in selected
            and candidate_extension(
                candidate.neighbors,
                source_edges,
                coloring[: len(base_points)],
                color_count,
            )
            is None
        ]
        scan_seconds = time.monotonic() - scan_started
        if not blockers:
            exhausted = True
            break

        candidate = blockers[0]
        translation = circle_translation(candidate.r, branch=1)
        shifted = translated_points(source_points, translation.point)
        structured_cross = enumerate_structured_cross_edges(
            base_points,
            source_points,
            candidate.r,
            shifted,
            translation,
        )
        expected_cross = {
            (base_index, source_index)
            for source_index, base_indices in candidate.neighbors.items()
            for base_index in base_indices
        }
        if structured_cross != expected_cross:
            raise AssertionError("candidate neighbor table and cross edges disagree")

        mapped: list[int] = []
        added_vertices = 0
        overlap_vertices = 0
        for point in shifted:
            vertex = point_indices.get(point)
            if vertex is None:
                vertex = len(points)
                point_indices[point] = vertex
                points.append(point)
                added_vertices += 1
            else:
                overlap_vertices += 1
            mapped.append(vertex)
        for left, right in source_edges:
            if mapped[left] == mapped[right]:
                raise AssertionError("a translated source edge collapsed to a loop")
            edges.add(
                (min(mapped[left], mapped[right]), max(mapped[left], mapped[right]))
            )
        for base_index, source_index in structured_cross:
            translated_index = mapped[source_index]
            if base_index == translated_index:
                raise AssertionError("a structured cross edge collapsed to a loop")
            edges.add(
                (
                    min(base_index, translated_index),
                    max(base_index, translated_index),
                )
            )

        selected.append(candidate.source_r_index)
        solve_started = time.monotonic()
        encoding = ColoringSAT(
            len(points),
            edges,
            color_count,
            break_color_symmetry=False,
        )
        next_coloring = encoding.solve()
        solve_seconds = time.monotonic() - solve_started
        iterations.append(
            {
                "iteration": iteration_index,
                "blocking_candidates_before_add": len(blockers),
                "selected_source_r_index_zero_based": candidate.source_r_index,
                "second_neighbor_count": candidate.second_neighbor_count,
                "vertices_added": added_vertices,
                "overlap_vertices_merged": overlap_vertices,
                "vertices_after_add": len(points),
                "listed_edges_after_add": len(edges),
                "scan_seconds": scan_seconds,
                "solve_seconds": solve_seconds,
                "five_sat_after_add": next_coloring is not None,
            }
        )
        if next_coloring is None:
            unsat = True
            break
        coloring = next_coloring

    if unsat:
        status = "PRIMARY_D4_LISTED_EDGE_GRAPH_5_UNSAT_WITHOUT_PROOF"
        coloring_payload: list[int] | None = None
    elif exhausted:
        status = "PRIMARY_D4_POOL_EXHAUSTED_WITH_5_COLORING"
        coloring_payload = list(coloring)
    else:
        status = "PRIMARY_D4_ITERATION_LIMIT_WITH_5_COLORING"
        coloring_payload = list(coloring)

    if coloring_payload is not None and not ColoringSAT(
        len(points),
        edges,
        color_count,
        break_color_symmetry=False,
    ).validate(coloring):
        raise AssertionError("final CEGIS coloring is invalid")

    _write_json(
        output_graph_path,
        serialize_graph(
            "2510-point closure plus primary D4 translated Heule copies",
            points,
            edges,
        ),
    )
    report: dict[str, object] = {
        "status": status,
        "lower_bound_6_proved": False,
        "candidate_pool": "primary-d4-conjugate-lengths-v1",
        "candidate_count": len(pool),
        "selected_source_r_indices_zero_based": selected,
        "iterations": iterations,
        "vertices": len(points),
        "listed_exact_unit_edges": len(edges),
        "whole_graph_5_sat": not unsat,
        "whole_graph_coloring": coloring_payload,
        "output_graph_path": str(output_graph_path),
        "elapsed_seconds": time.monotonic() - started,
        "edge_completeness_scope": (
            "complete within the original base, within every translated copy, "
            "and from the original base to each copy; extra unit edges between "
            "different translated copies are intentionally omitted"
        ),
    }
    _write_json(output_report_path, report)
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-graph",
        type=Path,
        default=Path(
            "artifacts/cegis/"
            "heule1061-second-closure-recertified-complete.graph.json"
        ),
    )
    parser.add_argument(
        "--base-report",
        type=Path,
        default=Path(
            "artifacts/cegis/heule1061-second-closure-recertified.json"
        ),
    )
    parser.add_argument(
        "--source-vtx",
        type=Path,
        default=Path("third_party/CNP-SAT/vtx/529.vtx"),
    )
    parser.add_argument(
        "--source-edge",
        type=Path,
        default=Path("third_party/CNP-SAT/edge/529.edge"),
    )
    parser.add_argument(
        "--output-graph",
        type=Path,
        default=Path(
            "artifacts/cegis/heule2510-plus-translated-heule529.graph.json"
        ),
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        default=Path(
            "artifacts/cegis/heule2510-plus-translated-heule529.json"
        ),
    )
    parser.add_argument("--r-vertex", type=int, default=98)
    parser.add_argument("--branch", type=int, choices=(0, 1), default=1)
    parser.add_argument(
        "--primary-d4-cegis",
        action="store_true",
        help="run the 24-member structured CEGIS pool instead of one copy",
    )
    parser.add_argument("--max-iterations", type=int, default=24)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.primary_d4_cegis:
        report = run_primary_d4_cegis(
            base_graph_path=args.base_graph,
            base_report_path=args.base_report,
            source_vertex_path=args.source_vtx,
            source_edge_path=args.source_edge,
            output_graph_path=args.output_graph,
            output_report_path=args.output_report,
            max_iterations=args.max_iterations,
        )
    else:
        report = run_single_copy(
            base_graph_path=args.base_graph,
            base_report_path=args.base_report,
            source_vertex_path=args.source_vtx,
            source_edge_path=args.source_edge,
            output_graph_path=args.output_graph,
            output_report_path=args.output_report,
            r_vertex=args.r_vertex,
            branch=args.branch,
        )
    summary = {
        key: value
        for key, value in report.items()
        if key not in {"whole_graph_coloring"}
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
