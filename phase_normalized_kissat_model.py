"""Verify a Kissat model after reversing a recorded CNF phase transform."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from coloring_sat import Edge, parse_dimacs_edge
from kissat_coloring_model import parse_kissat_output
from lean_coloring_cnf import decode_lean_coloring_model
from phase_normalized_cnf import unflip_model


def validate_phase_normalized_kissat_coloring(
    *,
    log_path: str | Path,
    mapping_path: str | Path,
    vertex_count: int,
    edges: Iterable[Edge],
    colors: int = 5,
) -> tuple[tuple[int, ...], dict[str, int | str]]:
    """Reverse the sign bijection, decode one color, and check every edge."""

    status, transformed_model = parse_kissat_output(log_path)
    if status != "SAT":
        raise ValueError(f"Kissat result is {status}, not SAT")
    mapping = json.loads(Path(mapping_path).read_text(encoding="utf-8"))
    raw_flipped = mapping.get("true_variables")
    if not isinstance(raw_flipped, list) or not all(
        isinstance(variable, int) and variable > 0
        for variable in raw_flipped
    ):
        raise ValueError("phase mapping has no valid true_variables list")
    original_model = unflip_model(transformed_model, raw_flipped)
    normalized_edges = {
        (min(left, right), max(left, right)) for left, right in edges
    }
    coloring = decode_lean_coloring_model(
        original_model,
        vertex_count=vertex_count,
        edges=normalized_edges,
        colors=colors,
    )
    return coloring, {
        "status": "VALIDATED_SAT_COLORING",
        "transformed_model_literals": len(transformed_model),
        "flipped_variables": len(set(raw_flipped)),
        "vertex_count": vertex_count,
        "edge_count": len(normalized_edges),
        "color_count": colors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--edge", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--colors", type=int, default=5)
    args = parser.parse_args()

    vertex_count, edges = parse_dimacs_edge(args.edge)
    coloring, stats = validate_phase_normalized_kissat_coloring(
        log_path=args.log,
        mapping_path=args.mapping,
        vertex_count=vertex_count,
        edges=edges,
        colors=args.colors,
    )
    payload: dict[str, object] = {
        **stats,
        "source_log": str(args.log),
        "phase_mapping": str(args.mapping),
        "model_literals": stats["transformed_model_literals"],
        "class_sizes": [
            coloring.count(color) for color in range(args.colors)
        ],
        "coloring": list(coloring),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(
        json.dumps(
            {key: value for key, value in payload.items() if key != "coloring"},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
