# Source provenance

Extracted from the Julia implementation on
`HomogeneousTools/PartialFlagVarieties.jl`, branch
`feat/atomic-decomposition`, commit
`eaad8dd2f339e29671692c0ef336e2ef46cc0b56`.
The input was the later exact-rational snapshot preserved in
`ZeroLociSearch.jl/reproducibility/quantum/PartialFlagVarieties/src/GromovWitten.jl`.
Its SHA-256 is
`05e3d8ce7b542568d6ed74bd4b7ead8123617ad172509631a0f6fbcf782b6c62`.

That snapshot applies `exact-rational-blocks.patch`: quantum blocks are summed
in `Rational{BigInt}` at two torus points. The new package preserves this
arithmetic and removes the unused modular scalars and reconstruction helpers.
No Python runtime or copied PartialFlagVarieties implementation is required.
The dependency is the registered PartialFlagVarieties 0.4 package.

The mathematical source identified by the Julia implementation is
Alexeev, Cavenaghi, Galindo, Georgiev, Katzarkov, and Martins,
*On the atomic decomposition of complete intersection in flag varieties*,
[arXiv:2609.21914v2](https://arxiv.org/abs/2609.21914v2).
The inherited fourfold fixtures retain their example numbers.
The free-tree and canonical coloring iterators retain the GKMtools attribution
and MIT notice in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

The extraction splits Schubert data, bundle weights, tree enumeration,
localization, pair-block evaluation, exact linear algebra, Novikov polynomials,
and quantum multiplication into separate source files. It also fixes automatic
enumeration in zero Chern directions, duplicate explicit curve classes,
degree-zero primary invariants with more than three marks, and validates
Novikov polynomial construction. These changes have dedicated regression tests.

A subsequent review replaced marked-tree Cartesian enumeration with independent
vertex sums, added cached Schubert subword automata and normalized Euler
factors, and reused exact matrix factorizations. It also fixed shared Schubert
word and free-tree iterator storage, explicit detection of rational poles,
and duplicate normalized Novikov monomials. Small subword and marked-orbit
enumerations check the optimized formulas independently; inherited matrices
remain regression fixtures. The internal manual describes the conventions
and ownership rules, and `benchmark/README.md` records measured comparisons.

Historical computation snapshots remain unchanged so that their recorded
hashes and reproducibility remain valid. To migrate existing scripts, replace
the patched dependency with registered PartialFlagVarieties and add:

```julia
using PartialFlagVarieties
using GromovWittenFlags
```

The computational API names are preserved. A package version containing the
old embedded quantum implementation should be replaced before loading both
packages, to avoid conflicting exports.
