# Fano fourfold comparison, 27 October 2026

Keep the current GromovWittenFlags computation as the production backend.
GKMtools is useful as an independent validation implementation and has broader
GKM functionality. The measured Euler-twisted fourfold workload does not justify
adding it and Oscar as required dependencies.

The distinction between its public matrix API and its integration kernel
matters: retaining GKMtools' classical pairing makes the five smaller examples
slightly faster than our public API, but its prepared computation is slower on
the nonlinear Küchle examples. These are workload-specific measurements.

## Method

Apple M4, 32 GiB RAM, macOS arm64, Julia 1.13.0, one Julia thread.
GromovWittenFlags 0.1.0 at commit `4a8f864`, PartialFlagVarieties 0.4.0,
Semisimple 1.3.0, Oscar 1.8.2, and GKMtools 1.0.0 at
[`a19904b`](https://github.com/mgemath/GKMtools.jl/tree/a19904bdb5b9c0158163b356dffb1cb79fd7a3c8).
Package loading is excluded; the recorded process used
`--compiled-modules=existing` after normal Oscar precompilation failed locally.

Each entry is a median of three serial, warmed samples. All classical metrics
and every individual Novikov coefficient matrix agreed exactly, after warm-up
and after every timed sample. The common basis is the nondegenerate ambient
quotient selected by GromovWittenFlags, transported to independently constructed
native GKMtools models. The default uses Billey Schubert lifts in both libraries.
Every operator is the flag-ambient projection of multiplication by the target's
first Chern class; additional zero-locus cohomology lies outside this comparison.
See [README.md](README.md) for the full method, adapter, and reproduction commands.

Raw measurements:
[matrix API](results-matrix.toml),
[prepared integration](results-prepared.toml),
[Küchle degree one](results-kuechle-degree1.toml),
[K3-29 by degree](results-k3-by-degree.toml),
[Schur-lift control](results-schur.toml), and
[Weyl checks](results-weyl.toml).
The full dependency resolution is in [Manifest-v1.13.toml](Manifest-v1.13.toml).

`GKMtools prepared` retains its classical matrix and inverse pairing. Its timed
positive-degree work uses the public unmarked integration routine with the Euler
twist and the same dimension filter as its matrix API. Its additional preparation
is excluded from the prepared timings and recorded separately. Our public API
still includes classical work and two exact torus evaluations. This comparison
therefore favors the prepared GKMtools path.

## Complete smaller operators

Times in milliseconds. The complete curve-degree lists are 0–1 for the cubic
and 0–2 for the other four examples. All targets have complex dimension four.

| Fourfold | Basis size | GromovWittenFlags | GKMtools matrix API | GKMtools prepared |
|:--|--:|--:|--:|--:|
| Cubic in P⁵ | 5 | 1.088 | 27.763 | 0.533 |
| Gushel–Mukai, O(1) ⊕ O(2) on Gr(2,5) | 6 | 7.237 | 146.712 | 6.441 |
| Quartic in P⁵ | 5 | 3.672 | 37.948 | 2.594 |
| Quadric–cubic in P⁶ | 5 | 6.934 | 102.752 | 5.228 |
| Three quadrics in P⁷ | 5 | 11.415 | 402.796 | 9.406 |

Prepared GKMtools is 1.1–2.0× faster for these repeated computations, with fewer
Julia-reported allocations. Its additional classical preparation took about
18–244 ms, depending on the example. Combined graph construction, basis transport,
and metric verification are excluded from every column. Those combined setup
times cannot be attributed to either library individually.

## Nonlinear bundles on Grassmannians

These are truncations through ambient curve degree three, including degree zero.
The full c5 and b8 operators extend through degree five; they were checked in our
extended regression suite, but are not full cross-library benchmarks here.

| Fourfold and defining bundle | GromovWittenFlags (s) | GKMtools prepared (s) | Additional GKMtools preparation (s) |
|:--|--:|--:|--:|
| c5: Λ²S* ⊕ Λ³Q ⊕ O(1) on Gr(3,7) | 2.326 | 27.065 | 19.324 |
| b8: O(1)³ ⊕ Sym²S* on Gr(2,7) | 0.767 | 6.018 | 3.454 |

Our API is about 11.6× and 7.8× faster, respectively, even with GKMtools' classical
preparation excluded. Julia-reported allocations were 2.00 GB versus 30.87 GB
for c5 and 0.78 GB versus 4.97 GB for b8. These are cumulative allocated bytes,
not peak memory; native-library allocations may be omitted.

A separate public-matrix comparison through degree **one** gave 0.0393 s versus
28.589 s for c5 and 0.0148 s versus 4.213 s for b8. It should not be confused with
the degree-three prepared comparison. Rebuilding the symbolic classical pairing
is a substantial part of that API's cost.

## Product fourfold K3-29

The model is the complete intersection O(1,2) ⊕ O(1,1) in P² × P⁴,
with first Chern degrees (1,2) and ambient quotient basis size nine.
The cross-library run covers all nine automatic curve classes of total
coordinate sum at most three, one class at a time. The full automatic list
contains twelve classes; (3,1), (4,0), and (5,0) were not included in this run.
The positive-degree measurements below are seconds:

| Curve class | GromovWittenFlags | GKMtools matrix API | GKMtools prepared |
|:--|--:|--:|--:|
| (0,1) | 0.00361 | 0.387 | 0.327 |
| (1,0) | 0.00279 | 0.237 | 0.160 |
| (0,2) | 0.00644 | 2.947 | 2.905 |
| (1,1) | 0.00617 | 1.380 | 1.324 |
| (2,0) | 0.00429 | 2.003 | 1.818 |
| (1,2) | 0.01728 | 18.046 | 17.510 |
| (2,1) | 0.01292 | 17.268 | 17.253 |
| (3,0) | 0.00834 | 19.297 | 18.691 |

Degree zero also agreed; the prepared path simply returns its retained matrix,
so its degree-zero timing measures a cache lookup. The native preparation took
about 0.238 s, in addition to combined model and transport setup.
These rows must not be added and presented as a single full-operator timing:
our API repeats its classical work for each row.

GKMtools' higher-rank edge-multiplicity routine constructs a module kernel and
a polyhedron, then enumerates interior lattice points. The rank-one path instead
uses weighted compositions. The initial attempt at a full K3-29 operator was
interrupted during its first GKMtools warm-up; its stack was inside the
polyhedral multiplicity routine. No full GKMtools median is available.
Our separate full-operator measurement was 0.201 s, covering all twelve classes.
The componentwise nonnegative degree bounds available in our flag model let
the specialized enumerator prune inadmissible assignments much earlier.
See the measured
[GKMtools multiplicity routine](https://github.com/mgemath/GKMtools.jl/blob/a19904bdb5b9c0158163b356dffb1cb79fd7a3c8/src/homology/curve_classes.jl).

## Equivariant-lift control

Using ordinary Schur polynomials as equivariant lifts gives the same ordinary
answers, but denser restrictions and substantially different performance.
For three quadrics, GKMtools' complete public matrix calculation took 4.295 s
with those lifts versus 0.403 s with its Billey lifts. We therefore used the
sparse Billey lifts for the main comparison and retained both raw datasets.
This control prevents the less favorable input representation from driving
the backend decision.

## Weyl-group overlap

PartialFlagVarieties uses Semisimple for its root systems, Weyl groups, and weight
lattices. GKMtools uses Oscar's corresponding objects in its generalized flag
implementation. Thus the overlap is chiefly **Semisimple versus Oscar**, rather
than another standalone Weyl implementation inside PartialFlagVarieties.

The representations differ. Semisimple encodes Dynkin type and rank in types,
with static root data; Oscar uses parent objects and its own algebraic element
types. There is no automatic shared Weyl-element type. The benchmark's explicit
reduced-word conversion passed for every element in A3, B2, C3, and G2:

| Type | Elements | Simple-root action checks | Mismatches |
|:--|--:|--:|--:|
| A3 | 24 | 72 | 0 |
| B2 | 8 | 16 | 0 |
| C3 | 48 | 144 | 0 |
| G2 | 12 | 24 | 0 |

Cartan matrices and Coxeter lengths agreed, and the transported words gave a
bijection in each case. This verifies these conventions and examples, rather
than asserting complete API equivalence or a performance ranking of Weyl groups.
See [weyl_checks.jl](weyl_checks.jl).

GKMtools' fast native type-A flag constructor stores block-permutation tuples;
it does not need Oscar Weyl elements to enumerate these Grassmannian vertices.
The benchmark uses that constructor. Its public Billey Schubert API instead
expects Weyl-backed vertices, which explains the small input adapter.
Replacing Semisimple throughout PartialFlagVarieties would also affect bundle
characters and highest-weight computations; it is a separate architectural
decision unsupported by these GW timings.

Source references at the measured revision:
[native type-A flags](https://github.com/mgemath/GKMtools.jl/blob/a19904bdb5b9c0158163b356dffb1cb79fd7a3c8/src/examples/homogeneous/standard.jl),
[generalized flags](https://github.com/mgemath/GKMtools.jl/blob/a19904bdb5b9c0158163b356dffb1cb79fd7a3c8/src/examples/homogeneous/GP.jl),
and [Semisimple's scope](https://github.com/HomogeneousTools/Semisimple.jl#relationship-to-oscar).

## Decision and limits

Keep GKMtools in this isolated benchmark environment. Retain the specialized
root-theoretic implementation for the current public package. Its advantages
on the nonlinear bundles and products, exact two-point checks, and existing
PartialFlagVarieties representation outweigh a modest repeated-integration win
on the five simpler examples. The prepared adapter remains useful for regression
comparisons and future experiments.

GKMtools' untwisted quantum Schubert engine reconstructs ordinary G/P quantum
products using different algorithms. It is not an Euler-twisted zero-locus
backend. Its generic GKM spaces, equivariant calculations, descendants, and
positive-genus features could motivate an optional integration for new scope;
they were not benchmarked here. See its
[quantum Schubert manual](https://github.com/mgemath/GKMtools.jl/blob/a19904bdb5b9c0158163b356dffb1cb79fd7a3c8/docs/src/GW/QuantumSchubert.md).

No claim is made about larger degrees, other geometries, multithreading,
startup on other Julia versions, or a future optimized GKMtools revision.
