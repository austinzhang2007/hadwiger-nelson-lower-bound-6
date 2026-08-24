"""Flip CNF variable signs so a supplied near-model is the all-false phase."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence


def canonicalize_residual_coloring(
    coloring: Sequence[int],
    *,
    first_residual_color: int = 3,
) -> tuple[int, ...]:
    """Swap the two residual colors if the higher one occurs first."""

    low = first_residual_color
    high = low + 1
    first_low = next(
        (index for index, color in enumerate(coloring) if color == low),
        len(coloring),
    )
    first_high = next(
        (index for index, color in enumerate(coloring) if color == high),
        len(coloring),
    )
    if first_high >= first_low:
        return tuple(int(color) for color in coloring)
    return tuple(
        high if color == low else low if color == high else int(color)
        for color in coloring
    )


def phase_true_variables(
    coloring: Sequence[int],
    *,
    colors: int = 5,
    first_residual_color: int = 3,
) -> set[int]:
    """Return true original and residual-prefix variables for one coloring."""

    canonical = canonicalize_residual_coloring(
        coloring, first_residual_color=first_residual_color
    )
    if any(color < 0 or color >= colors for color in canonical):
        raise ValueError("coloring contains a color outside the encoding")
    vertex_count = len(canonical)
    true_variables = {
        vertex * colors + color + 1
        for vertex, color in enumerate(canonical)
    }
    seen = False
    auxiliary_base = vertex_count * colors
    for vertex, color in enumerate(canonical):
        seen = seen or color == first_residual_color
        if seen:
            true_variables.add(auxiliary_base + vertex + 1)
    return true_variables


def rewrite_cnf_for_false_phase(
    source_path: str | Path,
    target_path: str | Path,
    true_variables: Iterable[int],
) -> dict[str, int]:
    """Apply the literal bijection x -> not x for the supplied variables."""

    source = Path(source_path)
    target = Path(target_path)
    flipped = {int(variable) for variable in true_variables}
    if any(variable < 1 for variable in flipped):
        raise ValueError("flipped variables must be positive")
    variables = declared_clauses = None
    parsed_clauses = 0
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    with source.open(encoding="ascii") as input_handle, temporary.open(
        "w", encoding="ascii"
    ) as output_handle:
        for line_number, raw_line in enumerate(input_handle, 1):
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("c"):
                output_handle.write(raw_line)
                continue
            fields = stripped.split()
            if fields[0] == "p":
                if len(fields) != 4 or fields[1] != "cnf":
                    raise ValueError(f"malformed header at line {line_number}")
                variables, declared_clauses = int(fields[2]), int(fields[3])
                if flipped and max(flipped) > variables:
                    raise ValueError("flipped variable exceeds CNF header")
                output_handle.write(raw_line)
                continue
            if variables is None or fields[-1] != "0":
                raise ValueError(f"malformed clause at line {line_number}")
            literals = [int(raw) for raw in fields[:-1]]
            transformed = [
                -literal if abs(literal) in flipped else literal
                for literal in literals
            ]
            output_handle.write(" ".join(map(str, transformed)))
            output_handle.write(" 0\n")
            parsed_clauses += 1
    if variables is None or declared_clauses != parsed_clauses:
        raise ValueError("CNF header clause count does not match input")
    temporary.replace(target)
    return {
        "variables": variables,
        "clauses": parsed_clauses,
        "flipped_variables": len(flipped),
    }


def unflip_model(
    model: Sequence[int], true_variables: Iterable[int]
) -> tuple[int, ...]:
    """Map a transformed DIMACS model back to the source variable signs."""

    flipped = {int(variable) for variable in true_variables}
    return tuple(
        -literal if abs(literal) in flipped else int(literal)
        for literal in model
    )


def _load_coloring(path: Path) -> list[int]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for key in ("whole_graph_coloring", "coloring", "best_coloring"):
        raw = payload.get(key)
        if isinstance(raw, list):
            return [int(color) for color in raw]
    raise ValueError("coloring report has no coloring")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cnf", type=Path, required=True)
    parser.add_argument("--coloring-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--colors", type=int, default=5)
    args = parser.parse_args()

    coloring = canonicalize_residual_coloring(
        _load_coloring(args.coloring_report)
    )
    true_variables = phase_true_variables(coloring, colors=args.colors)
    stats = rewrite_cnf_for_false_phase(
        args.cnf, args.output, true_variables
    )
    expected_variables = len(coloring) * (args.colors + 1)
    if stats["variables"] != expected_variables:
        raise ValueError(
            "source CNF is not the expected residual-prefix coloring encoding"
        )
    payload: dict[str, object] = {
        **stats,
        "source_cnf": str(args.cnf),
        "source_coloring_report": str(args.coloring_report),
        "vertex_count": len(coloring),
        "color_count": args.colors,
        "residual_colors": [3, 4],
        "canonical_coloring": list(coloring),
        "true_variables": sorted(true_variables),
    }
    if args.metadata is not None:
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.metadata.with_suffix(args.metadata.suffix + ".tmp")
        temporary.write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        temporary.replace(args.metadata)
    print(
        json.dumps(
            {
                key: value
                for key, value in payload.items()
                if key not in {"canonical_coloring", "true_variables"}
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
