using PartialFlagVarieties
using GromovWittenFlags

X = projective_space(2)
pt = point_class(X)
@assert gromov_witten_invariant(X, [pt, pt], [1]) == 1
@assert gromov_witten_invariant(X, fill(pt, 5), [2]) == 1
result = small_quantum_multiplication(X)
@assert evaluate.(result.matrix) == [0 0 3; 3 0 0; 0 3 0]
println("Quantum multiplication by c₁ on ℙ²:")
display(result.matrix)
