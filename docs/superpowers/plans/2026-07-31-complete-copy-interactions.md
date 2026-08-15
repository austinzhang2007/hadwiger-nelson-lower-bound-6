# Complete Copy Interactions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reconstruct the complete 8329-vertex unit-distance graph, export three exact SAT instances, and obtain independently checkable SAT or UNSAT evidence.

**Architecture:** A focused `copy_interactions.py` module proves and enumerates all copy-to-copy edges using quadratic-basis separation. `coloring_sat.py` gains a constant-size clique symmetry breaker. A separate verifier reconstructs the same graph without importing search code, while a proof runner treats timeout as `UNKNOWN`.

**Tech Stack:** Python 3.13, SymPy exact radicals, PySAT/CaDiCaL/Glucose, DRAT-trim, pytest.

## Global Constraints

- No floating-point comparison may decide a retained or omitted certificate edge.
- Only `stage3` contains the complete graph on all 8329 coordinates.
- UNSAT is accepted only with an independently verified DRAT/LRAT trace.
- SAT models must be checked against every edge of their exact graph.
- The workspace is not a Git repository, so commit steps are recorded as unavailable rather than fabricated.

---

### Task 1: Exact complete copy interactions

**Files:**
- Create: `copy_interactions.py`
- Test: `tests/test_copy_interactions.py`

**Interfaces:**
- Produces: `complete_copy_interactions(source_points, selected_r_indices) -> InteractionResult`
- `InteractionResult.edges` contains `(left_copy, left_vertex, right_copy, right_vertex)`.
- `InteractionResult.pair_certificates` records all 55 pairs.

- [ ] Write a failing test asserting 55 pair certificates, exactly three nonempty pairs with counts 648, 655, 651, and 1954 total edges.
- [ ] Run the test and confirm failure because `copy_interactions` is missing.
- [ ] Implement same-extension cross-coordinate grouping and mixed-extension product-coefficient exclusion.
- [ ] Directly verify every returned edge has exact squared distance one.
- [ ] Run the focused test and confirm it passes.

### Task 2: Constant-size SAT symmetry breaking

**Files:**
- Modify: `coloring_sat.py`
- Modify: `tests/test_coloring_sat.py`

**Interfaces:**
- Extend `ColoringSAT(..., symmetry_clique: Sequence[int] = ())`.
- Every pair in `symmetry_clique` must be an edge; vertex `i` receives color `i`.

- [ ] Write a failing test on a triangle plus isolated vertex.
- [ ] Confirm failure because the constructor lacks `symmetry_clique`.
- [ ] Add clique validation and one unit clause per clique vertex.
- [ ] Confirm the focused and existing coloring tests pass.

### Task 3: Exact graph/EDGE/CNF stage exporter

**Files:**
- Modify: `copy_interactions.py`
- Create: `tests/test_copy_interaction_export.py`

**Interfaces:**
- Produces `export_interaction_stages(...) -> list[StageArtifact]`.
- Stages have 52026, 52681, and 53332 edges.
- Each CNF fixes triangle vertices `(548, 1149, 668)` zero-based to colors `(0,1,2)`.

- [ ] Write a failing export test using a temporary directory.
- [ ] Confirm the exporter is missing.
- [ ] Implement atomic JSON writes and existing DIMACS/CNF writers.
- [ ] Validate stage edge nesting and exact counts.
- [ ] Run the focused test.

### Task 4: Independent complete-graph verifier

**Files:**
- Create: `copy_interaction_verifier.py`
- Create: `tests/test_copy_interaction_verifier.py`

**Interfaces:**
- Produces `verify_complete_interaction_artifact(...) -> dict[str, object]`.
- Must not import `translated_copy_cegis` or `copy_interactions`.

- [ ] Write a failing test against exported stage3.
- [ ] Confirm failure because the verifier is missing.
- [ ] Independently reconstruct all 55 pair classes and expected edges.
- [ ] Check coordinate identity, duplicate count, complete edge equality, and optional coloring.
- [ ] Run the focused test.

### Task 5: Proof-oriented solver runner

**Files:**
- Create: `proof_solver.py`
- Create: `tests/test_proof_solver.py`

**Interfaces:**
- `solve_cnf(cnf_path, result_path, proof_path, solver_name, phase_path=None)`.
- Result status is exactly `SAT`, `UNSAT_UNVERIFIED`, or `UNKNOWN_INTERRUPTED`.

- [ ] Write failing tests for a tiny SAT CNF and tiny UNSAT CNF.
- [ ] Confirm the runner is missing.
- [ ] Implement model/proof serialization and atomic result writes.
- [ ] Validate the tiny UNSAT proof with bundled `drat-trim`.
- [ ] Run focused tests.

### Task 6: Run all three stages

**Files:**
- Create under: `artifacts/interactions/`
- Modify: `research_log.md`
- Modify: `README.md`

- [ ] Export and hash all stage graph, edge, CNF, and report files.
- [ ] Solve stage1 with saved 5-coloring as phases and independently validate SAT.
- [ ] Run stage2 and stage3 with at least two solver configurations; preserve proof output for UNSAT.
- [ ] If UNSAT, validate with DRAT-trim; if interrupted, record `UNKNOWN` exactly.
- [ ] Run the independent stage3 geometry verifier.
- [ ] Run `.venv/bin/python -m pytest -q` and record the fresh result.
- [ ] Update documentation with exact statuses, runtimes, commands, and hashes.
