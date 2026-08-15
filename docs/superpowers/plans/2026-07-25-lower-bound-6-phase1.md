# Hadwiger–Nelson Lower Bound 6 Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development or superpowers:executing-plans to
> implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Reproduce a public 5-chromatic unit-distance graph and execute one
exactly certified coloring-blocking CEGIS iteration.

**Architecture:** Exact geometry, SAT, certificate verification, and search
are separate modules. Numeric filtering never decides certificate edges.

**Tech Stack:** Python 3.13, SymPy 1.13, SciPy 1.15, NetworkX 3.4,
python-sat 1.9.dev7, drat-trim.

## Global Constraints

- Exact algebraic coordinates are mandatory for retained candidates.
- Complete unit-distance adjacency is rebuilt from coordinates.
- A solver `UNSAT` result is not accepted without DRAT/LRAT validation.
- The Moser lattice/ring is not the complete candidate universe.

### Task 1: Exact geometry and graph I/O

**Files:** `exact_geometry.py`, `tests/test_exact_geometry.py`

- [ ] Run the Moser and author-coordinate tests and observe missing-module
  failure.
- [ ] Implement restricted expression parsing, exact squared distance,
  complete-edge rebuilding, Mathematica import, and equal-radius unit-circle
  intersections.
- [ ] Run `pytest tests/test_exact_geometry.py -q`; require all tests to pass.

### Task 2: SAT coloring

**Files:** `coloring_sat.py`, `tests/test_coloring_sat.py`

- [ ] Observe failure before implementation.
- [ ] Implement exactly-one coloring CNF, edge clauses, precedence symmetry
  breaking, incremental enumeration, augmentation, and DIMACS output.
- [ ] Require Moser 3-UNSAT/4-SAT and Heule 529 4-UNSAT/5-SAT.

### Task 3: Independent verification

**Files:** `graph_verifier.py`, `tests/test_graph_verifier.py`

- [ ] Observe failure before implementation.
- [ ] Implement a verifier that imports only `exact_geometry`, reconstructs
  every unit edge, compares complete edge sets, and checks vertex uniqueness.
- [ ] Verify author CNF/DRAT with the separately built `drat-trim`.

### Task 4: Candidate generation and CEGIS

**Files:** `cegis.py`, `tests/test_cegis.py`

- [ ] Observe failure before implementation.
- [ ] Generate numerical unit-circle intersections, deduplicate, retain
  degree-at-least-five candidates, recover exact coordinates, and certify
  every reported base neighbor.
- [ ] Enumerate canonical 5-colorings, compute candidate coverage, choose a
  best candidate/greedy group, augment the graph, and obtain the next SAT
  coloring or a proof-producing UNSAT result.
- [ ] Save a deterministic JSON report under `artifacts/cegis/`.

### Task 5: Final audit

**Files:** `README.md`, `research_log.md`

- [ ] Run the full test suite and CLI verification commands from a clean
  process.
- [ ] Record exact counts, hashes, timings, coverage, new coloring status,
  and every failed branch.
- [ ] State clearly whether the result is a numerical experiment, a finite
  certificate for a known graph, or a new lower-bound proof.

