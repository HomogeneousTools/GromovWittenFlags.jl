# GromovWittenFlags.jl

GromovWittenFlags.jl is a **largely automatic Julia port of
[gwflags](https://github.com/leonardofcavenaghi/atomicdecompositions)**, the
Python/SageMath software accompanying
[*On the atomic decomposition of complete intersection in flag varieties*](https://arxiv.org/abs/2609.21914v2)
by Boris Alexeev, Leonardo F. Cavenaghi, Giovane Galindo, Bogdan Georgiev,
Ludmil Katzarkov, and Pedro Antonio Muniz Martins.
The translation and much of the subsequent restructuring, testing, and
optimization were carried out automatically using language-model tools.
Upstream authorship and license notices are preserved in the package's
`PROVENANCE.md` and `THIRD_PARTY_NOTICES.md`.

The package computes genus-zero primary Gromov–Witten invariants of partial
flag varieties and Euler-twisted invariants of smooth homogeneous zero loci.
It depends on PartialFlagVarieties.jl for the geometry and on Semisimple.jl
for root systems and Weyl groups.

All invariants are exact rational numbers. Small quantum multiplication is
returned as a matrix of sparse [`NovikovPolynomial`](@ref)s.

```@meta
DocTestSetup = quote
  using PartialFlagVarieties, GromovWittenFlags
end
```

```jldoctest
julia> X = projective_space(2);

julia> pt = point_class(X);

julia> gromov_witten_invariant(X, [pt, pt], [1])
1//1
```

Start with [Getting started](@ref), then read [Conventions and scope](@ref)
before interpreting zero-locus computations.
