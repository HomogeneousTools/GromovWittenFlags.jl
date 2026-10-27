# Implementation

The package separates root and Schubert data, bundle weights, free trees and
their symmetries, scalar localization, shared pair blocks, and quantum algebra.
It imports public APIs from PartialFlagVarieties rather than copying geometry
or extending its module at runtime.

The localization formulas come from the root-theoretic construction in
[arXiv:2609.21914v2](https://arxiv.org/abs/2609.21914v2).
Schubert restrictions use a cached automaton for the reduced-subword formula.
Decorated trees enumerate fixed loci; automorphism and cover factors are
included. Quantum blocks sum labeled mark placements on unmarked trees using
orbit–stabilizer, sharing graph enumeration between insertion pairs and the
two exact torus evaluations. Sparse-support pairs use rooted trees instead.
Covered-edge and low-valence vertex factors are cached locally. Endpoint Euler
factors are combined at the vertices before multiplication, reducing the size
of intermediate exact fractions. Synchronous callbacks reuse graph buffers.
Balanced streaming sums limit the sizes of intermediate rational denominators
without changing the exact arithmetic.

For each finite degree list, exact linear algebra reduces the ambient metric
and raises an index of the divisor pairing using one reusable LU factorization.
The divisor axiom reduces
positive-degree first-Chern insertions to two-point blocks.
The [Internals](@ref) page derives these choices and documents every private
helper, including signs, symmetry factors, ownership, and pole handling.

The default tests check classical and enumerative values, complete matrices
for small targets, non-simply-laced roots, products, bundle twists, input
validation, exact algebra, and exhaustive small-tree coloring orbits.
The six smaller paper fourfold matrices and Küchle matrices through degree three run
by default. `GWFLAGS_EXTENDED_TESTS=true` adds the full degree-five Küchle
c5 and b8 fixtures. These are regression data; they complement the independent
small enumerative checks.

The upstream gwflags copyright and permission notice and the GKMtools MIT
notice are retained alongside this package's MIT license. `THIRD_PARTY_NOTICES.md`
identifies the adapted components, and `licenses/gwflags-LICENSE` preserves the
upstream license unchanged. `PROVENANCE.md` records the original authors,
translation history, Julia source commit, exact-arithmetic patch, and snapshot
checksum. Historical reproducibility snapshots are preserved.
