"""Exact sqrt(7) rotation experiment outside the primary D4 field.

All pairs are tested by integer coefficient arithmetic in
Q(sqrt(3),sqrt(5),sqrt(7),sqrt(11)); no floating-point edge prefilter is used.
"""
import argparse
import hashlib
import json
from math import lcm, prod
from pathlib import Path
import time

import sympy as sp
from pysat.solvers import Solver

from coloring_sat import ColoringSAT, write_dimacs_edge
from exact_geometry import Point, parse_mathematica_vertices, serialize_graph


PRIMES = (3, 5, 7, 11)
RADICANDS = tuple(prod(p for bit, p in enumerate(PRIMES) if mask & (1 << bit))
                  for mask in range(16))


def rotate_exact(point, c, s):
    if sp.simplify(c*c + s*s - 1) != 0:
        raise ValueError("rotation is not exactly unit")
    return Point(c * point.x - s * point.y, s * point.x + c * point.y)


def rotate_sqrt7(point):
    return rotate_exact(point, sp.Rational(1, 8), 3 * sp.sqrt(7) / 8)


def coefficients(expr, radicands=RADICANDS):
    symbols = sp.symbols("r0:16")
    expanded = sp.expand(sp.radsimp(expr))
    replaced = expanded.xreplace({sp.sqrt(d): symbols[i]
                                 for i, d in enumerate(radicands) if i})
    values = [replaced.subs({s: 0 for s in symbols})]
    values.extend(replaced.coeff(symbols[i]) for i in range(1, 16))
    if any(v.is_Rational is not True for v in values):
        raise ValueError("coordinate outside the declared multiquadratic field")
    reconstructed = values[0] + sum(values[i] * symbols[i] for i in range(1, 16))
    if sp.expand(replaced - reconstructed) != 0:
        raise ValueError("coordinate outside the declared multiquadratic field")
    return tuple(values)


def integer_complete_edges(points, progress=None, primes=PRIMES):
    if len(primes) != 4 or len(set(primes)) != 4 or not all(sp.isprime(p) for p in primes):
        raise ValueError("field basis requires four distinct primes")
    radicands = tuple(prod(p for bit, p in enumerate(primes) if mask & (1 << bit))
                      for mask in range(16))
    raw = [(coefficients(p.x, radicands), coefficients(p.y, radicands)) for p in points]
    denominator = 1
    for xy in raw:
        for vector in xy:
            for q in vector:
                denominator = lcm(denominator, int(q.q))
    coordinates = [tuple(tuple(int(q.p) * (denominator // int(q.q)) for q in vector)
                         for vector in xy) for xy in raw]
    if len(set(coordinates)) != len(coordinates):
        raise ValueError("exact duplicate vertices")
    target = denominator * denominator
    edges = set()
    for u, left in enumerate(coordinates):
        for v in range(u + 1, len(coordinates)):
            right = coordinates[v]
            square = [0] * 16
            for axis in (0, 1):
                delta = [(i, a - b) for i, (a, b) in enumerate(zip(left[axis], right[axis]))
                         if a != b]
                for k, (i, a) in enumerate(delta):
                    square[0] += a * a * radicands[i]
                    for j, b in delta[k + 1:]:
                        square[i ^ j] += 2 * a * b * radicands[i & j]
            if square[0] == target and not any(square[1:]):
                edges.add((u, v))
        if progress and u % 100 == 0:
            progress(u + 1, len(points), len(edges))
    return edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("third_party/CNP-SAT/vtx/529.vtx"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--rotation", choices=("sqrt7", "anchor5over3", "anchor4over3"), default="sqrt7")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    base = parse_mathematica_vertices(args.source)
    c, sine, primes = {
        "sqrt7": (sp.Rational(1, 8), 3*sp.sqrt(7)/8, (3, 5, 7, 11)),
        "anchor5over3": (sp.Rational(7, 10), sp.sqrt(51)/10, (3, 5, 11, 17)),
        "anchor4over3": (sp.Rational(5, 8), sp.sqrt(39)/8, (3, 5, 11, 13)),
    }[args.rotation]
    points = list(base)
    seen = set(points)
    for point in base:
        candidate = rotate_exact(point, c, sine)
        if candidate not in seen:
            points.append(candidate)
            seen.add(candidate)
    def progress(done, total, edges):
        print(json.dumps({"stage": "all_exact_pairs", "vertices_scanned": done,
                          "vertices": total, "edges": edges}), flush=True)
    edges = integer_complete_edges(points, progress, primes=primes)
    graph_path = args.output / "rotation.graph.json"
    graph_path.write_text(json.dumps(serialize_graph("Heule529 + " + args.rotation, points, edges),
                                    indent=2) + "\n")
    write_dimacs_edge(args.output / "rotation.edge", len(points), edges)
    n = len(base)
    base_edges = {(u, v) for u, v in edges if v < n}
    encoding = ColoringSAT(len(points), edges, 5, break_color_symmetry=False)
    coloring = encoding.solve()
    sample_count = blocked = 0
    with Solver(name="cadical195", bootstrap_with=encoding.clauses) as solver:
        for sample in ColoringSAT(n, base_edges, 5).enumerate_colorings(limit=args.samples):
            sample_count += 1
            if not solver.solve(assumptions=[encoding.variable(v, c)
                                             for v, c in enumerate(sample)]):
                blocked += 1
    result = {
        "status": "EXACT_COMPLETE_GRAPH_SAT" if coloring else "UNSAT_WITHOUT_PROOF",
        "source": str(args.source), "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "coordinate_field": "Q(" + ",".join(f"sqrt({p})" for p in primes) + ")",
        "rotation": {"name": args.rotation, "cos": str(c), "sin": str(sine)},
        "vertices": len(points), "edges": len(edges), "base_edges": len(base_edges),
        "cross_edges": sum(u < n <= v for u, v in edges),
        "all_pairs_exactly_checked": len(points) * (len(points) - 1) // 2,
        "base_samples": sample_count, "blocked_samples": blocked,
        "coloring": coloring, "coloring_valid": coloring is not None and encoding.validate(coloring),
        "class_sizes": [coloring.count(c) for c in range(5)] if coloring else None,
        "lower_bound_6_proved": False, "elapsed_seconds": time.monotonic() - start,
    }
    (args.output / "report.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "coloring"}), flush=True)


if __name__ == "__main__":
    main()
