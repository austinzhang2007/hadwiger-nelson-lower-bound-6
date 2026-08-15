# Degree-3 Configuration CEGIS Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox syntax for tracking.

**Goal:** Continue phase-1 CEGIS with all exactly recovered degree-at-least-three
unit-circle candidates and arbitrary interacting blocking configurations.

**Architecture:** A new module builds candidate interaction edges with numeric
prefiltering followed by exact certification, extracts deletion-minimal
list-coloring blockers, resumes the saved phase-1 state, and alternates blocker
addition with full graph-coloring SAT. Certificate verification remains
independent.

**Tech Stack:** Python 3.13, SymPy 1.13, SciPy 1.15, python-sat 1.9.dev7,
drat-trim.

## Global Constraints

- Numeric comparisons may only prefilter; every retained unit relation is exact.
- A whole-pool extension is a finite failure witness only after all graph edges
  and coloring conflicts are independently checked.
- Solver UNSAT is not a lower-bound certificate until its DRAT/LRAT verifies.
- The deterministic seed remains `20260725`.

---

### Task 1: Candidate interaction graph

**Files:**
- Create: `configuration_cegis.py`
- Create: `tests/test_configuration_cegis.py`

**Interfaces:**
- Consumes: `cegis.CertifiedCandidate`
- Produces: `certified_candidate_unit_edges(candidates, tolerance=1e-9)`

- [ ] Write a test with three candidates where two pairs are exactly unit and a
  nearby non-unit pair is rejected.
- [ ] Run
  `.venv/bin/python -m pytest tests/test_configuration_cegis.py -q` and observe
  import failure for the missing module.
- [ ] Implement KD-tree pair prefiltering and retain a pair only when
  `left.point.squared_distance(right.point) == 1`.
- [ ] Rerun the focused test and require it to pass.

### Task 2: General blocking-core extraction

**Files:**
- Modify: `configuration_cegis.py`
- Modify: `tests/test_configuration_cegis.py`

**Interfaces:**
- Produces:
  `deletion_minimal_blocker(candidate_neighbors, candidate_edges, coloring, k)`

- [ ] Write a failing test for the V-shaped three-point configuration with
  endpoint lists `{4}`, `{3}` and center list `{3,4}`.
- [ ] Implement full-pool extension followed by deletion minimization through
  `cegis.candidate_extension`.
- [ ] Assert the returned three vertices are blocking and every one-vertex
  deletion is extendable.

### Task 3: Resume phase-1 state and execute CEGIS

**Files:**
- Modify: `configuration_cegis.py`
- Modify: `tests/test_configuration_cegis.py`

**Interfaces:**
- Produces:
  `run_configuration_cegis(coordinate_path, edge_path, phase1_report_path,
  output_path, max_iterations=100, minimum_candidate_degree=3)`

- [ ] Write a failing exact-point mapping test that rejects a phase-1 candidate
  absent from the new pool.
- [ ] Implement report parsing, exact coordinate mapping, selected-set recovery,
  blocker addition, augmented graph SAT, whole-pool extension witness output,
  and deterministic JSON serialization.
- [ ] If augmented SAT returns UNSAT, write CNF and DRAT beside the JSON and run
  the configured independent checker before setting any certificate flag.
- [ ] Run the focused tests.

### Task 4: Full experiment and independent verification

**Files:**
- Create: `artifacts/cegis/heule529-degree3-configurations.json`
- Create when applicable:
  `artifacts/cegis/heule529-degree3-configurations-full-pool.graph.json`
- Modify: `README.md`
- Modify: `research_log.md`

- [ ] Run the second-stage CLI with seed `20260725` and up to 100 iterations.
- [ ] For SAT termination, reconstruct all unit edges from serialized coordinates
  and check the entire explicit coloring; for UNSAT termination, verify DRAT.
- [ ] Record candidate/edge counts, every blocker core, augmented graph sizes,
  termination reason, elapsed time and SHA-256 hashes.
- [ ] Run `.venv/bin/python -m pytest -q` and Python byte-compilation.
