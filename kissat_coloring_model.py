"""Parse a Kissat log and independently validate a lean coloring model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from coloring_sat import Edge, parse_dimacs_edge
from lean_coloring_cnf import decode_lean_coloring_model


def parse_kissat_output(path: str | Path) -> tuple[str, tuple[int, ...]]:
    """Return a terminal Kissat status and all DIMACS model literals."""

    status: str | None = None
    model: list[int] = []
    with Path(path).open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            fields = raw_line.split()
            if not fields:
                continue
            if fields[0] == "s":
                if len(fields) != 2:
                    raise ValueError(f"malformed status at line {line_number}")
                status = {
                    "SATISFIABLE": "SAT",
                    "UNSATISFIABLE": "UNSAT",
                    "UNKNOWN": "UNKNOWN",
                }.get(fields[1])
                if status is None:
                    raise ValueError(
                        f"unknown status at line {line_number}: {fields[1]}"
                    )
            elif fields[0] == "v":
                for raw_literal in fields[1:]:
                    literal = int(raw_literal)
                    if literal:
                        model.append(literal)
    if status is None:
        raise ValueError("Kissat log has no terminal status")
    if status == "SAT" and not model:
        raise ValueError("SAT Kissat log has no model")
    if status != "SAT" and model:
        raise ValueError("non-SAT Kissat log unexpectedly contains a model")
    return status, tuple(model)


def validate_kissat_lean_coloring(
    log_path: str | Path,
    *,
    vertex_count: int,
    edges: Iterable[Edge],
    colors: int = 5,
) -> tuple[int, ...]:
    """Decode one color per vertex and validate every supplied graph edge."""

    status, model = parse_kissat_output(log_path)
    if status != "SAT":
        raise ValueError(f"Kissat result is {status}, not SAT")
    return decode_lean_coloring_model(
        model,
        vertex_count=vertex_count,
        edges=edges,
        colors=colors,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--edge", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--colors", type=int, default=5)
    args = parser.parse_args()

    vertex_count, edges = parse_dimacs_edge(args.edge)
    coloring = validate_kissat_lean_coloring(
        args.log,
        vertex_count=vertex_count,
        edges=edges,
        colors=args.colors,
    )
    payload = {
        "status": "VALIDATED_SAT_COLORING",
        "source_log": str(args.log),
        "vertex_count": vertex_count,
        "edge_count": len(edges),
        "color_count": args.colors,
        "class_sizes": [
            sum(color == selected for color in coloring)
            for selected in range(args.colors)
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
