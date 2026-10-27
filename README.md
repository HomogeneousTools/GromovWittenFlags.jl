# GromovWittenFlags.jl

Exact genus-zero Gromov–Witten invariants of partial flag varieties and smooth
zero loci of supported homogeneous bundles, using root-theoretic localization.

Part of [HomogeneousTools](https://homogeneous.tools), built on
[PartialFlagVarieties.jl](https://github.com/HomogeneousTools/PartialFlagVarieties.jl)
and [Semisimple.jl](https://github.com/HomogeneousTools/Semisimple.jl).
Varieties and bundles belong to PartialFlagVarieties; this package owns the
Schubert localization and quantum computations.

## Installation

Requires Julia ≥ 1.11 and PartialFlagVarieties 0.4. This package is currently
unregistered. From a local checkout:

```julia
using Pkg
Pkg.develop(path="/path/to/GromovWittenFlags.jl")
```

## Quick start

```julia
using PartialFlagVarieties
using GromovWittenFlags

X = projective_space(2)
pt = point_class(X)
gromov_witten_invariant(X, [pt, pt], [1])  # 1//1: a line through two points
gromov_witten_invariant(X, fill(pt, 5), [2])  # 1//1: a conic through five points

S = zero_locus(line_bundle(projective_space(3), 3))
gromov_witten_invariant(S, [], [1])  # 27//1: lines on a cubic surface

result = small_quantum_multiplication(X)
evaluate.(result.matrix)  # [0 0 3; 3 0 0; 0 3 0], all Novikov variables = 1
```

## Scope and conventions

- Exact `Rational{BigInt}` answers, checked at two nonsingular torus points.
- Schubert classes are minimal Weyl coset representatives; their lengths are
  complex codimensions. Ambient curve coordinates follow `marked_nodes(X)`.
- Zero-locus invariants use restricted ambient insertions and Euler twisting.
  Smoothness and regularity of the defining section are assumptions.
  Every defining weight must pair nonnegatively with the nonparabolic positive
  coroots; zero weights require a regularization that is not implemented.
- Quantum multiplication is the flag-ambient projection of multiplication by
  the **zero locus's** first Chern class.
  The result contains `matrix`, `grading`, `basis`, and `metric`. For zero loci,
  the basis is reduced modulo the radical of the ambient intersection pairing;
  additional cohomology is outside this calculation. The ambient subspace
  need not be preserved by the full quantum product.
- Explicit `betas` specify a finite truncation. Automatic enumeration requires
  strictly positive first Chern coordinates, or returns the zero operator when
  all coordinates vanish. Mixed zero and positive coordinates require explicit
  curve classes.

See the [manual](docs/src/index.md), [conventions](docs/src/conventions.md),
[API reference](docs/src/api.md), and runnable [examples](examples).
The extraction history and preserved notices are in
[PROVENANCE.md](PROVENANCE.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Development

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
julia devtools/format.jl --check
```

The default suite includes independent enumerative values, full quantum
matrices for six smaller fourfold examples, Küchle matrices through degree three,
exhaustive small-tree symmetry checks, and Aqua package checks.
Run the full inherited fourfold regressions with:

```sh
GWFLAGS_EXTENDED_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```

This additionally checks the full degree-five Küchle c5 and b8 matrices and
can take substantially longer. The same full suite runs in CI when the test
workflow is dispatched manually; ordinary pushes and pull requests run the
default suite.
See [CONTRIBUTING.md](CONTRIBUTING.md) for structure and review expectations.
The [internal manual](docs/src/internals.md) explains the localization factors,
automorphisms, cache ownership, and private helpers. Reproducible performance
measurements are in [benchmark/README.md](benchmark/README.md).
The [GKMtools comparison](benchmark/gkmtools/RESULTS.md) records fourfold timings,
Weyl-group compatibility checks, and the dependency decision.
See [REVIEW.md](REVIEW.md) for the correctness fixes and validation evidence.

## License

MIT; see [LICENSE](LICENSE) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
