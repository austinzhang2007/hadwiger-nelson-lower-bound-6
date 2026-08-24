"""Scan an entire projective D4 checkpoint for exact forced candidate pairs."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
from typing import Callable, Iterable, Mapping, Sequence

from coloring_sat import ColoringSAT, parse_dimacs_edge
from d4_blocker_cegis import (
    _projective_from_payload,
    _projective_payload,
    _write_json,
)
from d4_exact import D4ProjectivePoint, primary_d4_coordinate_backend
from d4_pair_search import (
    ForcedPoint,
    ForcedPointStats,
    exact_forced_pair_edges,
    exact_four_color_forced_points,
)
from exact_geometry import parse_mathematica_vertices
from graph_verifier import load_graph_json


def build_pair_scan_payload(
    *,
    source_report: str | Path,
    missing_color: int,
    candidates: Sequence[ForcedPoint],
    pair_edges: Iterable[tuple[int, int]],
    stats: ForcedPointStats,
    timings_seconds: Mapping[str, float],
) -> dict[str, object]:
    """Build the sparse portable certificate consumed by the apply step."""

    normalized_edges = sorted(
        {tuple(sorted((int(left), int(right)))) for left, right in pair_edges}
    )
    used_indices = sorted({index for edge in normalized_edges for index in edge})
    records: list[dict[str, object]] = []
    for index in used_indices:
        candidate = candidates[index]
        records.append(
            {
                "candidate_index": index,
                "forced_color": candidate.forced_color,
                "neighbors_zero_based": list(candidate.neighbors),
                "generator_pairs_zero_based": [
                    list(pair) for pair in candidate.generator_pairs
                ],
                "projective_certificate": {
                    **_projective_payload(candidate.projective),
                    "fifth_neighbor_zero_based": -1,
                },
            }
        )
    return {
        "status": "EXACT_PROJECTIVE_FORCED_PAIR_SCAN",
        "source_report": str(source_report),
        "missing_color": int(missing_color),
        "forced_candidates": len(candidates),
        "pair_edges": [list(edge) for edge in normalized_edges],
        "used_candidates": records,
        "search_stats": asdict(stats),
        "timings_seconds": {
            str(key): float(value) for key, value in timings_seconds.items()
        },
    }


def scan_projective_forced_pairs(
    *,
    base_graph_path: str | Path,
    primary_report_path: str | Path,
    source_vertex_path: str | Path,
    resume_report_path: str | Path,
    output_path: str | Path,
    missing_color: int,
    progress: Callable[[str, dict[str, object]], None] | None = None,
) -> dict[str, object]:
    """Recover all four-color-forced points and certify their unit pairs."""

    started = time.monotonic()
    base_points, _ = load_graph_json(base_graph_path)
    with Path(primary_report_path).open(encoding="utf-8") as handle:
        primary = json.load(handle)
    selected = [
        int(index) for index in primary["selected_source_r_indices_zero_based"]
    ]
    source_points = parse_mathematica_vertices(source_vertex_path)
    original_base_count = len(base_points) - len(selected) * len(source_points)
    affine = primary_d4_coordinate_backend(
        base_points[:original_base_count], source_points, selected
    )

    with Path(resume_report_path).open(encoding="utf-8") as handle:
        resume = json.load(handle)
    vertex_count, current_edges = parse_dimacs_edge(resume["edge_path"])
    candidate_projective: list[D4ProjectivePoint] = [
        affine.algebra.normalize_projective(
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
        for record in resume["iterations"]
    ]
    candidate_points = [affine.materialize(point) for point in candidate_projective]
    if vertex_count != len(base_points) + len(candidate_points):
        raise ValueError("resume vertex and certificate counts differ")
    coloring = tuple(int(color) for color in resume["whole_graph_coloring"])
    if not ColoringSAT(
        vertex_count, current_edges, 5, break_color_symmetry=False
    ).validate(coloring):
        raise ValueError("resume coloring is invalid")
    exact = affine.as_projective_backend(candidate_projective)
    all_points = [*base_points, *candidate_points]
    if progress is not None:
        progress(
            "checkpoint_loaded",
            {"vertices": vertex_count, "edges": len(current_edges)},
        )

    recovery_started = time.monotonic()
    candidates, stats = exact_four_color_forced_points(
        all_points,
        coloring,
        exact,
        missing_color=missing_color,
    )
    recovery_seconds = time.monotonic() - recovery_started
    if progress is not None:
        progress(
            "forced_points_recovered",
            {
                "candidates": len(candidates),
                "seconds": recovery_seconds,
                **asdict(stats),
            },
        )

    pair_started = time.monotonic()
    pair_edges = exact_forced_pair_edges(candidates, exact)
    pair_seconds = time.monotonic() - pair_started
    payload = build_pair_scan_payload(
        source_report=resume_report_path,
        missing_color=missing_color,
        candidates=candidates,
        pair_edges=pair_edges,
        stats=stats,
        timings_seconds={
            "forced_recovery": recovery_seconds,
            "pair_filter": pair_seconds,
            "total": time.monotonic() - started,
        },
    )
    _write_json(Path(output_path), payload)
    if progress is not None:
        progress(
            "scan_written",
            {
                "pair_edges": len(pair_edges),
                "used_candidates": len(payload["used_candidates"]),
                "seconds": pair_seconds,
                "output_path": str(output_path),
            },
        )
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-graph", required=True)
    parser.add_argument("--primary-report", required=True)
    parser.add_argument("--source-vertices", required=True)
    parser.add_argument("--resume-report", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--missing-color", required=True, type=int, choices=range(5))
    args = parser.parse_args()

    def report(event: str, details: dict[str, object]) -> None:
        print(json.dumps({"event": event, **details}, sort_keys=True), flush=True)

    result = scan_projective_forced_pairs(
        base_graph_path=args.base_graph,
        primary_report_path=args.primary_report,
        source_vertex_path=args.source_vertices,
        resume_report_path=args.resume_report,
        output_path=args.output,
        missing_color=args.missing_color,
        progress=report,
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "forced_candidates": result["forced_candidates"],
                "pair_edges": len(result["pair_edges"]),
                "used_candidates": len(result["used_candidates"]),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    _main()
