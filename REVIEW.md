# Review of the extracted implementation

This review covers the standalone Julia package after extraction.
[PROVENANCE.md](PROVENANCE.md) identifies the mathematical source and preserved
snapshot. The changes retain exact rational arithmetic and the public
computational API.

## Correctness and failure handling

- Public Schubert constructors now copy Weyl words. A caller can no longer
  mutate a cached representative through a returned class. The first free-tree
  iterator result also has independent storage; collecting an iterator no
  longer corrupts that first result as its successor state advances.
- Rational poles are checked explicitly. Julia can represent `1//0`, so
  exception handling around division alone did not reliably detect poles.
  Covered edges are checked before support pruning, and final outputs must be
  finite. Only poles retry; inconsistent torus evaluations and unrelated
  exceptions propagate immediately.
- Graph enumeration uses synchronous callbacks. A pole or callback exception
  cannot leave a blocked Channel producer holding graph buffers. Borrowed
  storage and read-only structural arrays have explicit contracts.
- Novikov construction combines equivalent normalized exponent keys, removes
  cancellation zeros, and owns its dictionary. Infinite coefficients and
  invalid or oversized exponents are rejected. Exponent addition and public
  curve/dimension bounds use checked arithmetic.
- Bundle validation describes the actual nonnegative weight–coroot condition.
  Failure is not reported as a proof of nonconvexity. Zero torus weights are
  rejected with the missing fiber-equivariant regularization identified.

The corresponding regressions are in `test/schubert.jl`, `test/trees.jl`,
`test/novikov.jl`, `test/invariants.jl`, and `test/localization.jl`.

## Structure and performance

Root actions, reduced-word roots, inverse edge lookups, and root–coroot
pairings are cached independently of torus points. Typed function barriers
separate arithmetic from the heterogeneous variety cache. Per-evaluation
caches belong to one target, avoiding cross-bundle reuse.

Schubert restrictions use dense reduced-prefix automata. Cover-degree search
streams feasible assignments with suffix pruning. Orbit–stabilizer factors
labeled mark placements into vertex sums, and insertion pairs share graph
enumeration and restrictions. Endpoint Euler factors are combined before
rational multiplication; low-valence factors are memoized. Streaming balanced
sums reduce repeated gcd operations against large prefix denominators.

The classical metric and divisor pairing share a localization pass.
Incremental exact elimination selects the ambient quotient basis. One LU
factorization is reused across curve coefficients, followed by sparse Novikov
assembly. No explicit metric inverse or dense polynomial matrix product is
needed.

[benchmark/README.md](benchmark/README.md) records reproducible workloads and
measured timings and allocations. The [internal manual](docs/src/internals.md)
explains the formulas, signs, symmetry factors, ownership, and private APIs.
Documenter checks all docstrings, and the tests require documentation on
private functions, types, and constants. Inline comments explain the steps
whose mathematical or ownership invariants are easy to overlook.

## Validation and limits

At extraction-review commit `4a8f864`, the default suite passed 2,564 checks
with bounds checking, including Aqua.
It compares Schubert restrictions against literal subword enumeration and
marked sums against explicit labeled placements with marked stabilizers.
Graph colorings are checked against brute force on complete and irregular
graphs; known star symmetries test the large-tree fallback. Complete smaller
fourfold matrices and exact Küchle coefficients through degree three complement
independent enumerative values.

The extended suite at that commit passed all 2,567 checks, including the complete degree-five
Küchle c5 and b8 matrices, with bounds checking enabled (21m34s locally).

Strict documentation, doctests, executable examples, and pinned formatting
checks also pass. Local checks use Julia 1.13.0. The CI matrix is configured
for Julia 1.11, 1.12, and stable Julia on Linux, plus stable Julia on macOS;
those additional jobs have not been run locally.

High-degree graph enumeration remains combinatorial. Exact agreement at two
torus points checks consistency and does not prove the localization formula.
The defining section must satisfy the documented geometric assumptions.
Quantum matrices describe the ambient projection; they do not supply missing
primitive cohomology. These limitations are stated in the public manual.

## Julia idiom and documentation follow-up

The follow-up removes boxed captures in four production methods. Both degree
searches now use named recursive workers with explicit scratch arguments;
the degree-orbit key consumes the chosen permutation without capturing a
reassigned binding; the pair-block support probe has its own evaluation binding.
The enumeration workloads in `benchmark/benchmarks.jl` isolate these costs.

The large-tree encoder returns its key and symmetry order together, replacing
a shared mutable counter. Torus evaluations no longer retain unused tangent
Euler and input-weight fields. Restriction tables live in the numeric layer;
the invariant layer selects the pair-block optimization. Redundant coordinate
validation and duplicate root–coroot arithmetic were removed.

Public docstrings now contain executable examples and describe return fields,
coordinate conventions, and ownership. Private contracts distinguish full
Schubert indices, block positions, local roots, and cover degrees. The internal
manual maps the call layers and explains scratch restoration, cache effects,
and callback lifetimes. The free-tree documentation accounts for the owned
output copy, and coloring documentation states its normalized-rooting precondition.

Benchmark tools live in dedicated modules with CLI guards, explicit configuration,
specialized measurement callables, and separate model, runner, and Weyl-check files.
Checkpoint serialization failures preserve the last completed file. Their isolated
tests exercise configuration, checkpoint failure handling, word transport, and
native/prepared coefficient agreement on a cubic and K3-29.

The follow-up default suite passes 2,595 checks with bounds checking, including
Aqua and independent finite-curve enumeration checks. Strict Documenter checks,
public and private docstrings, doctests, and all benchmark reference comparisons
pass. The degree-five extended-suite result above belongs to the earlier
extraction-review commit.

## Internal contracts and iterator follow-up

Explicit zero-denominator and nonfinite-result guards now signal a private
`_LocalizationPole`. Only that exception retries a torus pair; ordinary
`DivideError`s propagate immediately. Regression checks exercise successful
retry, arithmetic-error propagation, disagreement, and exhaustion of all sixteen
pairs, including fallback from the eight small-coordinate pairs to the original
eight larger-coordinate pairs.

Coloring successors reset complete preorder subtrees in place using stored
parent indices. The output-copy contract is preserved. Brute-force comparisons
now cover rooted colorings as well as unrooted colorings, with a distinguished
root mark enforcing root-preserving isomorphisms. An isolated workload measures
32.5% fewer allocated bytes; the method and scope are recorded in
`benchmark/README.md`. Rational row selection also converts directly into its
owned scratch vector, avoiding an intermediate collection.

Private docstrings distinguish written buffers from return values and describe
the lifetimes of rooted shape and coloring data. Novikov documentation specifies
rank, exponent, and merge contracts. Küchle matrix regressions compare exact
polynomials rather than their rendered strings; all four rewritten fixtures
were checked coefficient for coefficient against their previous forms.

The default suite passes 2,867 checks with bounds checking, including Aqua.
Strict Documenter checks, doctests, pinned formatting, and warmed benchmark
reference comparisons pass. The complete degree-five sums were not rerun in
this follow-up; their earlier validation remains recorded above.

## Certified backend and public 0.1.0

The release adds branch recursion, full convex Euler indices, formal denominator
bounds, integer magnitude bounds, deterministic modular recovery, process-local
workers and atomic source-bound checkpoints. These are opt-in computational
tools; ordinary rational tree localization remains the default. The certified
API and checkpoint schema are experimental in 0.1.0.

The candidate passed 502,590 package checks on Aion with Julia 1.11.4, including
Aqua, parallel worker cleanup, corruption and interruption recovery, and
certificate identity rejection. Generic certified comparisons cover a
Grassmannian, a cubic fourfold and a nonsplit presentation of the projective
plane. Branch checks also compare individual torus values with the independent
tree enumerator, including non-simply-laced root systems and multiple marks.
Three additional complete Küchle c5/b8 branch-matrix comparisons passed.

Strict Documenter checks and doctests, the pinned formatting check, and all
41 checks of the portable Peskine verifier passed locally. The latter reads a
completed result and checks its polynomial algebra export; it does not repeat
localization. The earlier full Peskine regeneration and separate checkpoint
replay are recorded as a historical milestone in `benchmark/peskine`.
Those measurements belong to the preserved prototype source, not to the
subsequently formatted release tree.

Case-specific geometry and reference values live under `examples/peskine`.
The arithmetic engine contains no Peskine model switch. Its certificates do
not establish the geometric localization hypotheses, closure of an ambient
projection, monodromy, or the passage from spectral factors to Hodge atoms.
