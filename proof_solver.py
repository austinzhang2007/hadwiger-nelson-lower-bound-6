"""Run a DIMACS CNF in a separate proof-oriented, result-recording process."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
from typing import Iterable, Sequence

from pysat.formula import CNF
from pysat.solvers import Solver, pysolvers


def decode_graph_coloring_model(
    model: Sequence[int],
    *,
    vertex_count: int,
    color_count: int,
    edges: Iterable[tuple[int, int]],
) -> tuple[int, ...]:
    """Decode exactly one positive color literal and validate every edge."""

    positive = {literal for literal in model if literal > 0}
    coloring: list[int] = []
    for vertex in range(vertex_count):
        selected = [
            color
            for color in range(color_count)
            if vertex * color_count + color + 1 in positive
        ]
        if len(selected) != 1:
            raise ValueError(
                f"vertex {vertex} has {len(selected)} selected colors"
            )
        coloring.append(selected[0])
    if any(coloring[left] == coloring[right] for left, right in edges):
        raise ValueError("model contains a monochromatic graph edge")
    return tuple(coloring)


def export_validated_coloring(
    *,
    result_path: str | Path,
    coloring_path: str | Path,
    phase_path: str | Path,
    vertex_count: int,
    color_count: int,
    edges: Iterable[tuple[int, int]],
) -> dict[str, object]:
    """Validate a saved SAT model and export a compact coloring and phase file."""

    result_path = Path(result_path)
    coloring_path = Path(coloring_path)
    phase_path = Path(phase_path)
    with result_path.open(encoding="utf-8") as handle:
        result = json.load(handle)
    if result.get("status") != "SAT":
        raise ValueError("solver result is not SAT")
    model = result.get("model")
    if not isinstance(model, list) or not all(
        isinstance(literal, int) and literal != 0 for literal in model
    ):
        raise ValueError("SAT result has no valid integer model")
    coloring = decode_graph_coloring_model(
        model,
        vertex_count=vertex_count,
        color_count=color_count,
        edges=edges,
    )
    class_sizes = [
        sum(color == selected for color in coloring)
        for selected in range(color_count)
    ]
    payload: dict[str, object] = {
        "status": "VALIDATED_SAT_COLORING",
        "source_result_path": str(result_path),
        "solver": result.get("solver"),
        "vertex_count": vertex_count,
        "color_count": color_count,
        "class_sizes": class_sizes,
        "coloring": list(coloring),
    }
    _write_json(coloring_path, payload)
    phases = [
        vertex * color_count + color + 1
        if color == coloring[vertex]
        else -(vertex * color_count + color + 1)
        for vertex in range(vertex_count)
        for color in range(color_count)
    ]
    _write_json(
        phase_path,
        {
            "source_coloring_path": str(coloring_path),
            "vertex_count": vertex_count,
            "color_count": color_count,
            "phases": phases,
        },
    )
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def _load_phases(path: Path | None) -> list[int]:
    if path is None:
        return []
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    raw = payload.get("phases") if isinstance(payload, dict) else payload
    if not isinstance(raw, list) or not all(
        isinstance(literal, int) and literal != 0 for literal in raw
    ):
        raise ValueError("phase file must contain a list of nonzero literals")
    return [int(literal) for literal in raw]


def solve_cnf(
    *,
    cnf_path: str | Path,
    result_path: str | Path,
    proof_path: str | Path | None,
    solver_name: str,
    phase_path: str | Path | None = None,
) -> dict[str, object]:
    """Solve one CNF and atomically record SAT, UNSAT, or interruption."""

    cnf_path = Path(cnf_path)
    result_path = Path(result_path)
    proof = Path(proof_path) if proof_path is not None else None
    phases = _load_phases(Path(phase_path) if phase_path is not None else None)
    formula = CNF(from_file=str(cnf_path))
    started = time.monotonic()
    payload: dict[str, object] = {
        "status": "RUNNING",
        "solver": solver_name,
        "cnf_path": str(cnf_path),
        "cnf_sha256": _sha256(cnf_path),
        "variables": formula.nv,
        "clauses": len(formula.clauses),
        "proof_path": str(proof) if proof is not None else None,
        "phase_path": str(phase_path) if phase_path is not None else None,
    }
    _write_json(result_path, payload)

    try:
        with Solver(
            name=solver_name,
            bootstrap_with=formula.clauses,
            with_proof=proof is not None,
        ) as solver:
            if phases:
                solver.set_phases(phases)
            satisfiable = solver.solve()
            elapsed = time.monotonic() - started
            if satisfiable:
                model = [
                    literal
                    for literal in solver.get_model()
                    if abs(literal) <= formula.nv
                ]
                payload.update(
                    {
                        "status": "SAT",
                        "elapsed_seconds": elapsed,
                        "model": model,
                    }
                )
            else:
                if proof is None:
                    raise ValueError(
                        "UNSAT result requires a proof output path"
                    )
                trace = solver.get_proof()
                if trace is None:
                    raise RuntimeError(
                        f"{solver_name} returned UNSAT without a proof"
                    )
                proof.parent.mkdir(parents=True, exist_ok=True)
                temporary = proof.with_suffix(proof.suffix + ".tmp")
                with temporary.open("w", encoding="ascii") as handle:
                    for line in trace:
                        handle.write(line.rstrip())
                        handle.write("\n")
                temporary.replace(proof)
                payload.update(
                    {
                        "status": "UNSAT_UNVERIFIED",
                        "elapsed_seconds": elapsed,
                        "proof_sha256": _sha256(proof),
                        "model": None,
                    }
                )
    except KeyboardInterrupt:
        payload.update(
            {
                "status": "UNKNOWN_INTERRUPTED",
                "elapsed_seconds": time.monotonic() - started,
                "model": None,
            }
        )
    except pysolvers.error as exc:
        if "Caught keyboard interrupt" not in str(exc):
            raise
        payload.update(
            {
                "status": "UNKNOWN_INTERRUPTED",
                "elapsed_seconds": time.monotonic() - started,
                "model": None,
            }
        )
    _write_json(result_path, payload)
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cnf", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--proof", type=Path)
    parser.add_argument("--solver", default="glucose4")
    parser.add_argument("--phases", type=Path)
    args = parser.parse_args(argv)
    result = solve_cnf(
        cnf_path=args.cnf,
        result_path=args.result,
        proof_path=args.proof,
        solver_name=args.solver,
        phase_path=args.phases,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] == "SAT":
        return 10
    if result["status"] == "UNSAT_UNVERIFIED":
        return 20
    return 130


if __name__ == "__main__":
    raise SystemExit(main())
