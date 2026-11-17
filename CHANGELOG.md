# Changelog

## 0.1.0

First public release of the standalone HomogeneousTools package, a **largely
automatic Julia port of [gwflags](https://github.com/leonardofcavenaghi/atomicdecompositions)**.
The translation and much of the subsequent development used language-model tools;
upstream authorship and MIT notices are preserved.

- Exact primary genus-zero flag and Euler-twisted invariants, and small
  quantum multiplication in the ambient quotient basis.
- Rooted-branch summation as an alternative to decorated-tree enumeration.
- `ConvexEulerTwist` for an actual convex bundle represented by its full
  graded Euler index, with explicit geometric hypotheses.
- An experimental certified branch backend: formal denominator bounds,
  integer magnitude bounds, adaptive CRT, parallel workers and source-bound
  atomic checkpoints.
- Public exact flag integration and cyclic algebra analysis utilities.
- Executable Peskine example and numerical benchmark, separate from the core.
- Preserved upstream provenance, licenses and language-model-assistance notice.

The package is unregistered; install the `v0.1.0` Git tag by URL. Ordinary
rational arithmetic remains the default. Certified APIs and checkpoint schema
are experimental. Preserve the engine and environment with an archive;
prototype checkpoint identities change when the source or Julia runtime changes.

Primary spectral factors are algebraic output. Geometric Hodge-atom claims
require additional hypotheses and are not inferred by the package.
