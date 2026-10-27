using PartialFlagVarieties
using GromovWittenFlags

P3 = projective_space(3)
S = zero_locus(line_bundle(P3, 3))
@assert gromov_witten_invariant(S, [], [1]) == 27

P4 = projective_space(4)
Q = zero_locus(line_bundle(P4, 2))
result = small_quantum_multiplication(Q)
@assert evaluate.(result.matrix) == [0 0 6 0; 3 0 0 6; 0 3 0 0; 0 0 3 0]
println("Ambient quantum multiplication by c₁ on a quadric threefold:")
display(result.matrix)
