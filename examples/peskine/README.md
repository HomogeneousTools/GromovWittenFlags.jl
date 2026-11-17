# Peskine numerical example

This example supplies the geometry to the general package engine; the engine
contains no Peskine model switch. It is intended as an executable research
case and a reproducibility checkpoint for later optimizations.

From the package root, with its dependencies installed:

```sh
julia --project=. examples/peskine/compute.jl OUTPUT_DIRECTORY 4
julia --project=. examples/peskine/verify.jl OUTPUT_DIRECTORY fresh
```

The first command starts four local Julia worker processes. Run it inside a
scheduler allocation for the full computation. Use `0` workers for serial
execution. A third argument restricts maximum intrinsic degree, for example
`OUTPUT_DIRECTORY 0 0` computes only classical initialization.
No old numerical file or saved basis is read by the computation.
The separate verifier reads the reference fixture only after the run finishes.

The model is the regular zero locus of an actual globally generated rank-21
bundle, together with two hyperplanes, on `Fl(1,4;10)`. Its equivariant graded
class represents the full index. Intrinsic degree `d` maps to `(d,3d)`.
The example rebuilds its six-dimensional ambient quotient basis and computes
the complete divisor operator through intrinsic degree five.

Outputs include the coefficient blocks, all required two-point invariants,
the metric and basis conversion, cyclic algebra and projectors, the native
universal-line trace and surface moments, and the primitive scalar assignment.
The verifier exports all six polynomial multiplication operators and checks
their algebra identities, including the reconstructed higher-degree products.

The actual bundle, convexity, regularity, curve-degree embedding and smooth
surface-of-lines model are geometric hypotheses. The rank-22 inherited and
rank-32 vanishing Hodge representations, absolute irreducibility and regular
F-manifold genericity argument are additional inputs to the coarse atom statement.
The program does not infer these hypotheses or construct associated surfaces.

Re-running the same model and source with the same output directory reuses
completed per-degree checkpoints. To preserve a completed timing run, copy its
checkpoint directory into a different output directory before replaying it.
Changing the engine source or Julia version rejects the old identity; use a
fresh directory or the original archived source.

See [the benchmark milestone](../../benchmark/peskine/README.md) for the measured
prototype computation and the reference policy. It is not part of ordinary CI.
