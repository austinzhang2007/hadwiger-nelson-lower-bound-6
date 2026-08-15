# Second Unit-Circle Closure CEGIS Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox syntax for tracking.

**Goal:** Run arbitrary-configuration CEGIS over the second unit-circle
intersection closure of the exact 1061-vertex graph.

**Architecture:** `closure_cegis.py` loads an exact graph and a saved coloring,
uses the tested generic configuration step, checkpoints every iteration, and
emits either an explicit whole-pool coloring or a verified UNSAT certificate.

**Tech Stack:** Python 3.13, SymPy 1.13, SciPy 1.15, python-sat 1.9.dev7,
drat-trim.

## Global Constraints

- All retained coordinates and unit edges are symbolically certified.
- Numeric KD-tree operations only prefilter.
- Checkpoints contain the current base coloring and selected pool indices.
- Lower bound 6 remains false without complete exact geometry and verified DRAT.

---

### Task 1: Exact base-state loader

**Files:**
- Create: `closure_cegis.py`
- Create: `tests/test_closure_cegis.py`

- [ ] Write a failing test that loads a small exact graph and coloring from JSON,
  validates them, and rejects a same-color declared edge.
- [ ] Implement `load_exact_colored_graph(graph_path, report_path)`.
- [ ] Run the focused test and require it to pass.

### Task 2: Checkpointed closure loop

**Files:**
- Modify: `closure_cegis.py`
- Modify: `tests/test_closure_cegis.py`

- [ ] Write a failing serialization test for a completed synthetic iteration.
- [ ] Implement candidate generation, exact candidate interactions, generic
  configuration steps, atomic per-iteration JSON checkpoints, and deterministic
  resume fields.
- [ ] On pool extension, serialize all listed exact edges and the full coloring.
- [ ] On UNSAT, rebuild complete exact edges, generate CNF/DRAT and verify it.

### Task 3: Real experiment and audit

**Files:**
- Create: `artifacts/cegis/heule1061-second-closure.json`
- Create: `artifacts/cegis/heule1061-second-closure-listed.graph.json`
- Modify: `README.md`
- Modify: `research_log.md`

- [ ] Run up to 100 CEGIS iterations from the verified 1061-point coloring.
- [ ] Attempt independent complete-edge reconstruction on the terminal graph;
  report separately if resource limits prevent it.
- [ ] Check every listed edge exactly and every edge of the explicit coloring.
- [ ] Run the entire test suite and byte compilation.
