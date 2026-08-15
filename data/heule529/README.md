# Heule 529 derived symmetry-broken CNF

The author repository commit
`bb414955a6ef5f49f7df2b245b1e778aa67c068a` contains
`proof/529-4-sbp.drat`, but omits the matching `cnf/529-4-sbp.cnf`.
Running that proof against `cnf/529-4.cnf` fails at proof line 30758.

Every other graph in the same repository uses the same three
color-symmetry-breaking unit clauses:

```text
1 0
6 0
23 0
```

`529-4-sbp.cnf` is mechanically derived from the author's `529-4.cnf` by
inserting exactly those clauses after the header and increasing the declared
clause count from 11209 to 11212:

```bash
awk 'NR==1 {
  print "p cnf 2116 11212"
  print "1 0"
  print "6 0"
  print "23 0"
  next
} {print}' third_party/CNP-SAT/cnf/529-4.cnf \
  > data/heule529/529-4-sbp.cnf
```

SHA-256:
`fbe86fcd2c57e984b3b56e653ed0ea5d85396cb03ad63311d447514c568f5dee`.
With this input, unmodified `drat-trim` reports `s VERIFIED`.

