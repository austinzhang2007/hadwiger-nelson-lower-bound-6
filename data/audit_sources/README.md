# Arithmetic provenance

`primary_field_audit_20260906_initial.py` is the exact arithmetic module loaded
by the first persisted 52,034-vertex GMP audit. SHA-256:
`0cd99bb69c9bc36bf378fc0df101b9046339103c499099bab32e3a3d1411f611`.

The current root-level module has the same arithmetic kernel and additionally
rejects decimal/scientific coefficient strings and floating-point inputs.
The initial module accepted those strings as exact GMP rationals; it performed
no floating-point edge computation. The snapshot is retained for reproducibility,
not as the preferred input parser for new certificates.

This source snapshot by itself certifies no graph. Consult the persisted audit
JSON for its actual status, input hashes, checked-edge count and scope.
