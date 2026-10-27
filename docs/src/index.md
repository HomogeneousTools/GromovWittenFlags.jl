# GromovWittenFlags.jl

GromovWittenFlags computes genus-zero primary Gromov–Witten invariants of
partial flag varieties and Euler-twisted invariants of smooth homogeneous
zero loci. It depends on PartialFlagVarieties.jl for the geometry and on
Semisimple.jl for root systems and Weyl groups.

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
