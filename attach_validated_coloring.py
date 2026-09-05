"""Attach an independently validated 5-coloring to a CEGIS graph report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coloring_sat import parse_dimacs_edge


def _write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    temporary.replace(path)


def attach_validated_coloring(
    *,
    target_report_path: str | Path,
    coloring_result_path: str | Path,
) -> dict[str, object]:
    """Recheck every listed edge, then atomically update the target report."""

    target_path = Path(target_report_path)
    result_path = Path(coloring_result_path)
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if result.get("status") != "VALIDATED_SAT_COLORING":
        raise ValueError("coloring result is not independently validated SAT")
    raw_coloring = result.get("coloring")
    if not isinstance(raw_coloring, list):
        raise ValueError("coloring result has no coloring")
    coloring = tuple(int(color) for color in raw_coloring)

    with target_path.open(encoding="utf-8") as handle:
        target = json.load(handle)
    vertex_count, edges = parse_dimacs_edge(str(target["edge_path"]))
    if vertex_count != int(target["vertices"]):
        raise ValueError("target vertex count disagrees with edge file")
    if len(edges) != int(target["edges"]):
        raise ValueError("target edge count disagrees with edge file")
    if len(coloring) != vertex_count or any(
        color < 0 or color >= 5 for color in coloring
    ):
        raise ValueError("coloring has invalid length or color")
    if any(coloring[left] == coloring[right] for left, right in edges):
        raise ValueError("coloring violates the target edge file")

    target["status"] = "EXACT_INTERACTIONS_AUGMENTED_WITH_5_COLORING"
    target["lower_bound_6_proved"] = False
    target["whole_graph_5_sat"] = True
    target["whole_graph_coloring"] = list(coloring)
    metadata_keys = (
        "solver",
        "source_log",
        "phase_mapping",
        "vertex_count",
        "edge_count",
        "color_count",
        "model_literals",
        "class_sizes",
    )
    target["validated_5_coloring"] = {
        "result_path": str(result_path),
        **{key: result[key] for key in metadata_keys if key in result},
    }
    _write_json(target_path, target)
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-report", type=Path, required=True)
    parser.add_argument("--coloring-result", type=Path, required=True)
    args = parser.parse_args()
    payload = attach_validated_coloring(
        target_report_path=args.target_report,
        coloring_result_path=args.coloring_result,
    )
    print(
        json.dumps(
            {
                "status": payload["status"],
                "vertices": payload["vertices"],
                "edges": payload["edges"],
                "whole_graph_5_sat": payload["whole_graph_5_sat"],
                "class_sizes": payload["validated_5_coloring"].get(
                    "class_sizes"
                ),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
