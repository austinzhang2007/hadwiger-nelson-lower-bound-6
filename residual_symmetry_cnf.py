"""Add a linear CNF symmetry break for the two colors left after a fixed K3."""

from __future__ import annotations

from pathlib import Path


def add_residual_color_precedence(
    source: str | Path,
    target: str | Path,
    *,
    vertex_count: int,
    color_count: int = 5,
    first_residual_color: int = 3,
) -> dict[str, int]:
    """Require color 3 to occur before color 4 using prefix auxiliaries."""

    source = Path(source)
    target = Path(target)
    lines = source.read_text(encoding="ascii").splitlines()
    header_index = next(
        index for index, line in enumerate(lines) if line.startswith("p cnf ")
    )
    _, _, variables_raw, clauses_raw = lines[header_index].split()
    variables, clauses = int(variables_raw), int(clauses_raw)
    if variables != vertex_count * color_count:
        raise ValueError("CNF variable count is not the expected coloring encoding")
    color3, color4 = first_residual_color, first_residual_color + 1
    if color4 >= color_count:
        raise ValueError("two residual colors are required")
    added: list[str] = []
    for vertex in range(vertex_count):
        seen = variables + vertex + 1
        x3 = vertex * color_count + color3 + 1
        x4 = vertex * color_count + color4 + 1
        if vertex == 0:
            added.extend((f"-{x3} {seen} 0", f"-{seen} {x3} 0", f"-{x4} 0"))
            continue
        previous = variables + vertex
        added.extend(
            (
                f"-{previous} {seen} 0",
                f"-{x3} {seen} 0",
                f"-{seen} {previous} {x3} 0",
                f"-{x4} {previous} 0",
            )
        )
    lines[header_index] = (
        f"p cnf {variables + vertex_count} {clauses + len(added)}"
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text("\n".join([*lines, *added]) + "\n", encoding="ascii")
    temporary.replace(target)
    return {
        "original_variables": variables,
        "original_clauses": clauses,
        "variables": variables + vertex_count,
        "clauses": clauses + len(added),
        "added_variables": vertex_count,
        "added_clauses": len(added),
    }
