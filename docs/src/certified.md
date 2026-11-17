# Certified computation

The ordinary backend computes exact rational localization sums directly.
For larger degree rectangles, the branch backend replaces tree enumeration
with a triangular recurrence. Certified arithmetic evaluates that recurrence
over prime fields and proves uniqueness of the recovered rational answer.
All production arithmetic is Julia with the package's normal dependencies;
Python, SymPy and Oscar are not required.

## Selecting a backend

`algorithm=:branches` uses exact rational branch summation.
`arithmetic=:certified` always uses the branch recurrence, with proved bounds
and independent modular jobs. The result retains `matrix`, `grading`, `basis`
and `metric`, and additionally contains `certification` records.

```@meta
DocTestSetup = quote
    using PartialFlagVarieties, GromovWittenFlags
end
```

```jldoctest
julia> X = projective_space(1);

julia> result = small_quantum_multiplication(X; arithmetic=:certified);

julia> evaluate.(result.matrix) == [0 2; 2 0]
true

julia> all(c -> c["status"] == "certified", result.certification)
true
```

!!! note "Experimental in version 0.1.0"
    The opt-in certified API and checkpoint format may change in subsequent
    minor releases. Source-bound checkpoints deliberately reject an engine or
    runtime upgrade. Keep the original source and environment with a scientific
    archive, or regenerate its checkpoints after upgrading.

## Arithmetic guarantees and hypotheses

The current certified boundary is a dimension-zero **two-point block of
integral Schubert insertions**. Divisor equations convert those blocks to small
quantum multiplication. [`certified_pair_invariants`](@ref) exposes the same
block routine directly. Other numbers of marked points continue to use the
ordinary invariant API.

The denominator certificate uses Gauss valuations of formal localization
characters, before numerical torus specialization. The sharper valuation
recursion checks unit vertex-Euler contents and integral inverse-weight
contents. Other supported inputs use the coarse graph-factor bound.
This is a bound on the final constant rational integral, not on the numerical
denominators of individual torus contributions.

A second evaluation uses fixed-size, upward-rounded integer magnitude bounds.
After multiplying the invariant by its proved denominator, CRT recovery must
identify a unique integer in the proved range. No guessed denominator,
integrality assumption, floating-point rounding or probabilistic reconstruction
is used. Two nonsingular torus specializations must independently recover the
same rational values. As with ordinary localization, these arithmetic checks
do not establish the geometric localization or quantum-Lefschetz hypotheses.

The two default seeds are `(0,3)`. If exact preflight finds a singular seed,
supply another pair of distinct seeds in `0:31`. A pole appearing only modulo
a prime rejects that prime job and requests more good primes. Unrelated errors
are not classified as poles. `max_primes` bounds the prime search.

## Bundles represented by their full index

[`ConvexEulerTwist`](@ref) represents an **actual convex bundle** by the
equivariant K-class of its associated graded pieces. Its nonempty `certificate`
explains the actual bundle and convexity; a smooth regular zero locus is an
additional geometric hypothesis. Neither assertion is proved by the package.
Negative graded splittings contribute the full `H^0-H^1` edge index, with equal
formal characters cancelled before torus specialization. Surviving negative
multiplicities and identically zero characters are rejected.

This support is general: it is useful for filtered homogeneous bundles beyond
the Peskine model. The ordinary split-bundle validation remains unchanged.

## Workers and checkpoints

```julia
result = small_quantum_multiplication(X;
    arithmetic=:certified,
    workers=4,
    worker_heap_size="4G",
    checkpoint="checkpoints")
```

`workers=0` runs serially. A positive integer starts that many local Distributed
processes and cleans them up on exit; an existing vector of worker IDs uses
caller-owned processes. An active project containing GromovWittenFlags is
required. Submit the Julia driver within a suitable HPC allocation; the
package does not choose a scheduler, account, partition or memory request.
Each process owns its numerical caches and receives primitive job data rather
than live FLINT objects. One worker handles the magnitude bound, then the pool
handles independent prime/torus jobs. Curve degrees are currently sequential.
Internal fixed-point threading and sharing root coefficients across target
degrees are future optimizations.

Use one coordinator per checkpoint directory. It writes checksummed TOML
records using same-filesystem atomic replacement, saving each completed prime
job immediately. Per-degree directories contain the model manifest, the
denominator and magnitude record, individual prime jobs, and the final recovery.
Interrupted work can resume with a different worker count. Corrupt or
incompatible source, runtime, main arithmetic/geometry version, model,
insertion, curve-degree or torus identities are rejected. Checkpoints are
trusted local computation records, not adversarial proof files.

The versioned environment manifest records all registered dependencies.
The checkpoint identity also records the main arithmetic/geometry versions;
it is not a substitute for preserving the full manifest and source.
`worker_heap_size` is a per-process garbage-collection hint, not a hard limit.

## Classical integration and algebra analysis

[`equivariant_flag_integral`](@ref) integrates a homogeneous top-degree
characteristic expression or an array of such expressions at two torus points.
The callback receives actual bundle fiber characters and Schubert restrictions.
It is responsible for the degree condition. Supplied bundles may have zero
fiber weights because their Euler classes are not automatically divided out.

[`cyclic_quantum_algebra`](@ref) reconstructs a finite algebra from an exact
multiplication matrix and a cyclic unit. It computes characteristic/minimal
polynomials, all basis multiplication operators and rational primary projectors,
checking the algebra identities and an optional Frobenius metric.

To interpret this as a quantum algebra, the supplied operator must act on an
actual closed subalgebra. An ambient projection is not automatically closed.
Likewise, a small-point primary decomposition is not automatically a generic
Hodge-atom decomposition. Monodromy, Hodge equivariance and the relevant
deformation argument are case-specific mathematical inputs.

## Peskine as a reproducibility checkpoint

The executable example under `examples/peskine` supplies its flag/bundle model,
curve-degree embedding and geometric assertions. Those assumptions are not
hardcoded in the package. Its large measured computation and a separate
cross-process replay validate the numerical pipeline. The small reference
coefficients and benchmark summary are retained outside `src/` so future
optimizations can be compared with a concrete completed calculation.
