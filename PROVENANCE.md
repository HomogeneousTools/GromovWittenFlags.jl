# Source provenance

## Upstream gwflags

GromovWittenFlags.jl is a **largely automatic Julia port of
[gwflags](https://github.com/leonardofcavenaghi/atomicdecompositions)**, the
Python/SageMath implementation accompanying *On the atomic decomposition of
complete intersection in flag varieties* by Boris Alexeev, Leonardo F. Cavenaghi,
Giovane Galindo, Bogdan Georgiev, Ludmil Katzarkov, and Pedro Antonio Muniz Martins.
The translation and much of the subsequent restructuring, testing, and
optimization were carried out automatically using language-model tools.
The original gwflags implementation and accompanying mathematical work are
the upstream sources of this port.

The upstream [project credits](https://github.com/leonardofcavenaghi/atomicdecompositions/blob/affbf8ff6c251dc72a0bfaab5f3ddb35251e6d24/docs/about.md)
identify Ludmil Katzarkov and Maxim Kontsevich as project leaders and list
Alexeev, Cavenaghi, Galindo, Georgiev, and Martins as contributors.
Its MIT copyright notice names Leonardo Cavenaghi.
The [upstream license](https://github.com/leonardofcavenaghi/atomicdecompositions/blob/affbf8ff6c251dc72a0bfaab5f3ddb35251e6d24/LICENSE)
was checked at commit `affbf8ff6c251dc72a0bfaab5f3ddb35251e6d24` and is preserved
unchanged in [licenses/gwflags-LICENSE](licenses/gwflags-LICENSE).
The upstream file ends mid-sentence in its warranty disclaimer; the archived
copy intentionally preserves that text exactly.
This revision records the license check, not the Python revision used for the
initial translation. The package retains its own MIT license and the separate
GKMtools notice in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Julia extraction and development

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

## Certified branch extension

The 0.1.0 release incorporates the subsequent branch recursion, convex Euler
index twist, formal Gauss denominator certificates, integer magnitude bounds,
and modular recovery developed during the Peskine computation. Parallel job
orchestration and source-bound atomic checkpoints use Julia standard libraries.
The ordinary exact-rational API remains available and is the default.

These modules are general computational tools. The Peskine flag/bundle model,
its Hodge-theoretic hypotheses, reference numbers and performance milestone
are in `examples/peskine` and `benchmark/peskine`, outside production code.
The frozen prototype arithmetic and checkpoints remain in the separate original
computation archive. Formatting and release documentation change source-bound
checkpoint identities, so those prototype checkpoints must not be presented as
checkpoints produced by the final release source.

The extension was developed largely automatically using language-model tools,
like the earlier translation.
Its regression checks and exact numerical certificates complement review of
the mathematical formulas; they do not prove the geometric hypotheses.
