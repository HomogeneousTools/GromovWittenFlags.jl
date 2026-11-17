# GromovWittenFlags.jl

GromovWittenFlags.jl is a **largely automatic Julia port of
[gwflags](https://github.com/leonardofcavenaghi/atomicdecompositions)**, the
Python/SageMath software accompanying
[*On the atomic decomposition of complete intersection in flag varieties*](https://arxiv.org/abs/2609.21914v2)
by Boris Alexeev, Leonardo F. Cavenaghi, Giovane Galindo, Bogdan Georgiev,
Ludmil Katzarkov, and Pedro Antonio Muniz Martins.
The translation and much of the subsequent restructuring, testing, and
optimization were carried out automatically using language-model tools.
The original gwflags implementation and mathematical work are credited above;
their provenance and license notices are preserved below.

The package computes exact genus-zero Gromov–Witten invariants of partial flag
varieties and smooth zero loci of supported homogeneous bundles, using
root-theoretic localization.

Part of [HomogeneousTools](https://homogeneous.tools), built on
[PartialFlagVarieties.jl](https://github.com/HomogeneousTools/PartialFlagVarieties.jl)
and [Semisimple.jl](https://github.com/HomogeneousTools/Semisimple.jl).
Varieties and bundles belong to PartialFlagVarieties; this package owns the
Schubert localization and quantum computations.

## Installation

Requires Julia ≥ 1.11 and PartialFlagVarieties 0.4. This package is currently
unregistered. Install the public release directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/HomogeneousTools/GromovWittenFlags.jl", rev="v0.1.0")
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
  Ordinary bundles must have nonnegative weight–coroot pairings. An actual
  convex bundle with nonsplit graded pieces can instead be represented by
  `ConvexEulerTwist(graded; certificate="...")`; its geometric convexity and
  regular section are caller-supplied hypotheses. Zero weights require a
  regularization that is not implemented.
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

## Branch recursion and certified arithmetic

Use `algorithm=:branches` to resum decorated trees by a triangular power-series
recursion. The ordinary exact-rational API remains the default.
The opt-in `arithmetic=:certified` backend instead combines checked prime-field
evaluations with proved denominator and integer magnitude bounds. It returns
exact rational answers without assuming integrality. For example:

```julia
result = small_quantum_multiplication(projective_space(2);
    arithmetic=:certified, workers=2, checkpoint="gw-checkpoints")
```

The certified backend and checkpoint schema are **experimental in 0.1.0**.
The small test cases cover ordinary flags, split bundles and an index twist;
the complete Peskine computation is a larger regression benchmark. Supported
geometric hypotheses, singular-seed handling, checkpoint identity, and limits
are documented in the [certified computation manual](docs/src/certified.md).

`equivariant_flag_integral` provides exact classical flag integrals.
`cyclic_quantum_algebra` reconstructs a finite algebra when its unit is cyclic
for the supplied multiplication operator. Primary factors alone are not a
geometric atom decomposition; closure and Hodge-theoretic hypotheses must be
established separately.

The portable Peskine example is in [examples/peskine](examples/peskine).
Case-specific geometric assertions, numerical references and measured timings
remain outside the production modules. It is a reproducibility checkpoint for
future optimizations, not a special case in the quantum engine.

## License

MIT; see [LICENSE](LICENSE) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
The upstream gwflags copyright and permission notice is preserved in
[licenses/gwflags-LICENSE](licenses/gwflags-LICENSE).
