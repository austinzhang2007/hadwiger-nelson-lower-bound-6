# 2026-08-24—25 continuation log

This file records the post-move continuation while the desktop workspace root
still resolves through a symbolic link.  It is intended to be merged into
`research_log.md` once the canonical workspace binding permits patching
existing files.

Strict status labels remain unchanged:

- heuristic conflict minimization is a numerical/combinatorial experiment;
- a SAT model becomes a finite counterexample only after every listed edge is
  checked;
- no lower bound 6 claim is allowed without exact geometry and independently
  verified 5-UNSAT DRAT/LRAT.

## 49,820 vertices / 462,566 exact-listed edges

The persisted edge-ladder checkpoint had 438,605 active edges after 3,015 SAT
calls.  TabuCol seed `20260842`, sample size 8, ran 20,000,000 iterations in
472.46 seconds and reduced the target-edge conflict count from 95 to 16.
An independent standard-library recount parsed 49,820 colors and all 462,566
DIMACS edges and found exactly these 16 monochromatic edges (one-based):

```text
(223,9045) (246,25833) (301,809) (517,11790)
(632,9001) (1304,11589) (1311,8744) (1381,37997)
(1898,13292) (2460,32819) (8443,39076) (8780,21960)
(8828,34675) (9936,19769) (17483,19016) (17556,19111)
```

Three reproducible stagnation-kick Tabu trajectories gave best conflict counts
15, 16, and 16:

| seed | sample | kick | iterations | best |
|---:|---:|---:|---:|---:|
| 20260843 | 8 | 4 | 10,000,000 | 15 |
| 20260844 | 16 | 8 | 10,000,000 | 16 |
| 20260845 | 32 | 16 | 10,000,000 | 16 |

Frozen-neighbourhood exact SAT repair was locally UNSAT at radii 0, 1, 2,
and 3, with respectively 32, 1,123, 12,313, and 39,759 active vertices.
These results say only that the outside coloring cannot remain fixed.
They are not whole-graph UNSAT results.

An assumption-core repair run completed 50 fast rounds, relaxing 11,730
vertices, then reached a harder core.  Longer CaDiCaL195 and Glucose4 runs
were interrupted after about 30 and 25 minutes when a verified whole-graph
SAT model superseded them.  Their strict status is `UNKNOWN_INTERRUPTED`.

## Lean CNF and phase normalization

For graph coloring, the pairwise at-most-one clauses are unnecessary for
equisatisfiability: every vertex needs at least one true color, and adjacent
vertices may share no true color.  From any satisfying multi-color model,
choosing one true color per vertex gives an ordinary proper coloring.
`lean_coloring_cnf.py` and its tests implement this encoding.

- standard exact-one CNF: 249,100 variables, 2,860,853 clauses;
- lean CNF: 249,100 variables, 2,362,653 clauses;
- lean CNF plus residual 3/4 precedence: 298,920 variables,
  2,561,932 clauses.

The best-16 coloring was globally permuted so the fixed triangle
`(548,1149,668)` has colors `(0,1,2)`, and so color 3 occurs before color
4.  `phase_normalized_cnf.py` then flipped 99,637 variable signs, including
the residual prefix auxiliaries, so that this near-coloring is the all-false
phase.  Direct clause evaluation found exactly 16 unsatisfied clauses.

Kissat 4.0.4 was run with:

```bash
kissat --sat --phase=false --walkinitially=true --walkeffort=1000000 \
  --seed=2026082407 \
  heule49820-d4-blockers-5-lean-residual-sbp-tabu16-normalized-falsephase.cnf
```

It returned `s SATISFIABLE` / exit 10 after 2.64 seconds and 2,003
conflicts.  The independent verifier:

1. parsed all 298,920 model literals;
2. reversed the recorded sign bijection;
3. selected one true color per vertex from the lean model;
4. checked all 462,566 exact-listed edges.

The verified color-class sizes are
`[12354,10192,10556,11387,5331]`, with zero monochromatic listed edges.
The canonical compact coloring-array SHA-256 is
`869fb5968d676f24c71fc56677cd5ce856a3ab21cf87e5145546955d0e3d9e60`.
Therefore this graph is strictly **5-colorable** and cannot prove the lower
bound 6.

Artifact SHA-256 values:

| artifact | SHA-256 |
|---|---|
| target edge file | `9b95e64e20d331ae6a438acb8f4ce2bf5049bade12973a456d0a4f91db540877` |
| lean CNF | `120abe9aeeb59d07c90b6b59c817c8de55577f29c4bb19219cb64623a759def1` |
| lean residual CNF | `ac649ddf60240ef58851c85f6cb403f63acc067dd23ccd6171f8996276da82a2` |
| normalized false-phase CNF | `71ad9bb07ba0b2bb37f5dfe96e144c4cbb1a853a710da39e673e8e744dc2c1e1` |
| phase mapping | `f7e9c0aab91073019886633f2b0d9202a8d7ce8502f407d4abb12e7697b9d881` |
| Kissat log/model | `99031e122eb8be533af15a01c5c4af4d697548321fc804c33c4fb830e73af3b3` |
| independently verified coloring JSON | `1f3aed162bad8a67c57070c9cf35f6077ce531d9544b398e430c93f16a90e9c4` |

The 1.8 GB CEGIS report was atomically updated to
`EXACT_INTERACTIONS_AUGMENTED_WITH_5_COLORING`.

## Next CEGIS round

An all-projective forced-pair scan for missing color 4 was started from the
new coloring:

- base graph: `artifacts/cegis/heule2510-primary-d4-cegis.graph.json`;
- primary report: `artifacts/cegis/heule2510-primary-d4-cegis.json`;
- author source coordinates: `third_party/CNP-SAT/vtx/529.vtx`;
- resume report:
  `artifacts/cegis_49820_interactions/heule49820-d4-blockers.json`;
- seed dependence: none in exact recovery/pair certification.

The scan checked 99,921,538 cross-color center pairs and recovered 3,561
exact forced-color-4 candidates in 2,461.25 seconds.  Exact pair filtering took
241.87 seconds and certified 617 unit candidate pairs involving 977 candidates.
Applying all used candidates added 3,908 generator edges and all 617 pair edges,
producing a 50,797-vertex / 467,091-edge exact-listed graph.  CaDiCaL195 found
a 5-coloring in 4.40 seconds; an independent DIMACS-edge pass verified zero
conflicts and class sizes [12466,10289,10657,11496,5889].  The canonical
coloring-array SHA-256 is
`4a2da4445f27f504ab8f5abee96d4ae289b84f5e96aa06f53d15e610089daf76`.
A memory-bounded interaction-edge augmentation of this graph is now running.

## Current strict conclusion

- 50,797 / 467,091 exact-listed graph: **SAT**, independently validated;
- its 49,820 / 462,566 parent checkpoint is also **SAT**;
- complete induced unit-distance edge set: not yet certified for this large
  graph (the numerical annulus was only a prefilter);
- lower bound 6: **not proved**;
- 5-UNSAT DRAT/LRAT: none.
