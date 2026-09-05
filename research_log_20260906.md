# 2026-09-06 — field audit, independent arithmetic, and new rotations

Project conclusion: **lower bound 6 remains unproved**. There is no verified
5-UNSAT DRAT/LRAT certificate. This entry records actual finite experiments,
including negative results, and the precise scope of each audit.

## Recovered 52,034-vertex checkpoint

The previous missing-color-4 scan completed with 102,507,124 cross-color center
pairs, 4,076 exact forced points, and 831 exact unit candidate pairs involving
1,237 new vertices. Applying them added 4,948 generator edges and 831 internal
edges. CaDiCaL195 returned a five-coloring in 6.0366 seconds for the resulting
52,034-vertex / 478,193-listed-edge graph.

On resumption, an independent streaming read of the saved coloring and DIMACS
edge file checked all 478,193 edges: zero monochromatic edges, color classes
`[12540,10536,11012,11599,6347]`. Coloring-array SHA-256 (compact JSON) is
`b20eac6263ec8ad6393a1d8b7607ac9f3abdbd8954f3c265335765c116d0f1a9`.

Input checkpoint:
`artifacts/cegis_50797_next_projective_pairs_color4_applied/heule52034-d4-blockers.json`.
Its redundant radical graph JSON is about 1.44 GB; the report is about 2.59 GB.
The earlier full independent verifier runs left no successful persisted result
and no live process. Their geometry status is **unconfirmed**, not passed.

## Audit implementation and remaining gates

Read-only review identified that `d4_listed_graph_verifier.py` shares arithmetic
with the search, omits base edges from its unit checks, uses a floating-point
prefilter for duplicate detection, has no progress persistence, and builds a
large SAT CNF just to validate a saved coloring. It also reconstructs translated
base coordinates without comparing them against every redundant base JSON
coordinate. CPU activity alone did not locate its bottleneck precisely.

`primary_field_audit.py` now independently reduces the five defining radical
relations and performs rational arithmetic with GMP. Tests compare all 1,024
basis products with the existing implementation, four dense rational products
(seed 20260906), and exact homogeneous unit identities. Zero denominators and
near-unit nonedges are rejected. Coefficient parsing accepts only explicit
integers or fractions, rejecting decimal/scientific strings and floats.

`primary_checkpoint_audit.py` streams certificates with ijson and tests **all**
listed edges, including base edges. Every 10,000 checked edges it atomically
writes an audit JSON with input hashes and checked counts. It checks the saved
coloring directly without constructing CNF. The coordinate definition is the
original base prefix plus the author translation recipe plus homogeneous
certificates. It explicitly does **not** certify redundant large graph-JSON
consistency, duplicate absence, induced-edge completeness, or any UNSAT claim.
These remain final-certificate gates.

Current run (started this entry) writes:
`artifacts/cegis_50797_next_projective_pairs_color4_applied/heule52034-primary-audit.json`.
At this entry's checkpoint it has checked 370,000 edges with no invalid edges;
the run is still in progress. Read the persisted `status` and `checked_edges`
before interpreting later results.

The pre-existing exact interaction augmentation is also running from the same
52,034-point checkpoint, with output `artifacts/cegis_52034_interactions` and
`--skip-sat`. Its numerical annulus is a candidate prefilter, not a proof of
edge completeness. No new edge count is asserted until it writes its output.

## The field restriction and a primary-source check

The four-center recovery solves two perpendicular-bisector linear equations.
It therefore never leaves

\[
K=\mathbb Q(\sqrt3,\sqrt5,\sqrt{11},h_-,h_+),\qquad
h_\pm^2=(25\pm3\sqrt{33})/8.
\]

This is a commutative number field. “D4” refers to its nonabelian Galois
extension component, not to noncommutative coefficient multiplication.

The 32-term basis is independent: let B=Q(sqrt(3),sqrt(5),sqrt(11)). Neither
h-square is a square in B, because its norm under the involution negating
sqrt(3) is 41/8, which is not a square in the fixed biquadratic field
Q(sqrt(5),sqrt(11)). Their product 41/8 is also not a square in B (its square
class is 82). Thus the two square classes are independent and [K:B]=4.
This also justifies using a nonzero denominator coefficient vector as a
nonzero field element in the fixed primary field.

[Madore, arXiv:1509.07023, Proposition 4.6](https://arxiv.org/html/1509.07023)
gives a five-color upper bound for Q(sqrt(3),sqrt(11))^2 by reduction to F_11.
It does not directly give a five-coloring of K^2: in characteristic 11 the
new h-square reduces to 25/8 = 10, which is nonsquare in F_11. Consequently
the F_11 reduction used there does not extend as an F_11-valued reduction of
these integral generators. This is a limitation of that argument, **not** a
proof that K^2 requires six colors, nor a proof that no other upper-bound
method applies.

The [Moser lattice paper, arXiv:2606.12325](https://arxiv.org/abs/2606.12325)
likewise must not be applied to larger coordinate fields without verifying
its hypotheses.

To leave K explicitly, use new square roots at primes outside its ramification
set. The quartic polynomial for h is `8*t^4 - 50*t^2 + 41`; exact computation
gives discriminant `7406733312 = 2^11 * 3^6 * 11^2 * 41`. Together with the
base quadratic fields, K can ramify only at 2,3,5,11,41. Quadratic fields with
sqrt(7), sqrt(13), or sqrt(17) therefore cannot be subfields of K.

## Three completed domain-expansion experiments

`quadratic_rotation_experiment.py` takes the author's 529 coordinates, rotates
a copy about the origin, merges exact duplicates, and enumerates every pair
using integer coefficients in a four-prime multiquadratic field. All 558,096
pairs are checked; there is no floating-point edge prefilter in this experiment.
It solves the full graph and tests extension of 128 canonical base five-colorings.
Canonical means restricted-growth color precedence, removing color permutations;
the enumeration is deterministic and is not a uniform sample.

| rotation (cos, sin) | new prime | vertices | complete edges | extra edges beyond two copies | blocked / 128 |
|---|---:|---:|---:|---:|---:|
| (1/8, 3 sqrt(7)/8) | 7 | 1057 | 5340 | 0 | 0 |
| (7/10, sqrt(51)/10) | 17 | 1057 | 5376 | 36 | 0 |
| (5/8, sqrt(39)/8) | 13 | 1057 | 5361 | 21 | 0 |

All three complete graphs are **SAT**. Their color-class sizes are respectively
`[216,198,202,214,227]`, `[206,198,205,226,222]`, and
`[206,195,210,218,228]`. Runtime was approximately 27 seconds per experiment.
The first angle produces only two copies sharing the origin. The other two
were selected from author-vertex squared radii 5/3 and 4/3. The formula
cos(theta)=1-1/(2r^2) guarantees a unit chord from that vertex to its rotated
image and creates the extra certified interactions recorded above.

The `cross_edges` field in the machine report includes 36 origin-to-copy edges;
those already belong to a rotated copy. The table instead counts only edges
in excess of twice the original 2,670-edge graph, to avoid overstating interaction.

Published small exact graph, DIMACS edge, and coloring files are under
`data/rotation_experiments/`. Author coordinate SHA-256:
`ce0cf260e431972c1222521f1b552bc94e7c42040ddf4b4110cee0eaab518dbb`.
These results exclude the three particular augmented graphs as six-chromatic
witnesses. Zero coverage on 128 models is not a universal statement about all
base five-colorings or all possible rotated-copy configurations.

Reproduce, from the repository root:

```bash
.venv/bin/python quadratic_rotation_experiment.py --rotation sqrt7 --output artifacts/sqrt7_rotation_20260906 --samples 128
.venv/bin/python quadratic_rotation_experiment.py --rotation anchor5over3 --output artifacts/sqrt17_anchor_rotation_20260906 --samples 128
.venv/bin/python quadratic_rotation_experiment.py --rotation anchor4over3 --output artifacts/sqrt13_anchor_rotation_20260906 --samples 128
```

Verification this turn: 26 focused and related regression tests passed in
29.05 seconds. The test scope covers the new arithmetic, input validation,
all-edge audit persistence, rotation identities, complete integer edge scan,
and existing D4/coloring routines. This is not a full-suite rerun and does not
establish the research goal.

Next research step: collect the large audit and interaction outputs, solve the
augmented edge graph, and develop interacting configurations in the newly
introduced fields. A single independently rotated copy has not supplied a
blocker in the tested models.
