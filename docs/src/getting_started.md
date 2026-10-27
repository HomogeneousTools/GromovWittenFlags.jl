# Getting started

## Installation

For this unregistered package, develop a local checkout:

```julia
using Pkg
Pkg.develop(path="/path/to/GromovWittenFlags.jl")
```

Use `using PartialFlagVarieties, GromovWittenFlags` in computation scripts.
No patched branch of PartialFlagVarieties or Python environment is required.

## Invariants and zero loci

Schubert basis elements come from the ambient flag variety. The following
computes the number of lines on a smooth cubic surface:

```@meta
DocTestSetup = quote
  using PartialFlagVarieties, GromovWittenFlags
end
```

```jldoctest
julia> X = projective_space(3);

julia> S = zero_locus(line_bundle(X, 3));

julia> gromov_witten_invariant(S, [], [1])
27//1

julia> expected_gromov_witten_degree(S, [1], 0)
0
```

The alternative form is
`gromov_witten_invariant(X, [], [1]; bundle=line_bundle(X, 3))`.

## Quantum multiplication

```jldoctest
julia> X = projective_space(2);

julia> result = small_quantum_multiplication(X);

julia> evaluate.(result.matrix) == [0 0 3; 3 0 0; 0 3 0]
true

julia> evaluate(result.matrix[1, 3], [2])
6//1

julia> fano_curve_classes(X)
(3, [[0], [1]])
```

Columns are input basis elements and rows are output coefficients.
The `basis` field lists the corresponding ambient Schubert representatives;
`metric` is their intersection matrix and `grading` has diagonal entries
`codimension - dimension/2`.

Pass `betas=[[0], [1]]` to compute a specified finite set of curve classes.
Include the zero class to retain the classical term. Repeated classes are
deduplicated. Omitted classes contribute nothing to this truncation.
