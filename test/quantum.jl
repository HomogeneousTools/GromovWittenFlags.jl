@testset "Finite curve search versus bounded lattice enumeration" begin
  for coefficients in ([1, 2], [2, 3]), bound in 0:6
    actual = GWF._compact_beta_vectors(coefficients, bound)
    expected = [
      collect(beta) for beta in Iterators.product(0:bound, 0:bound)
      if sum(coefficients .* beta) <= bound
    ]
    sort!(expected; by=beta -> (sum(beta), beta))
    @test actual == expected
    @test length(unique(objectid.(actual))) == length(actual)
  end
end

@testset "Small quantum multiplication" begin
  for n in 1:3
    X = projective_space(n)
    result = small_quantum_multiplication(X)
    expected = zeros(QQ, n + 1, n + 1)
    for i in 1:n
      expected[i + 1, i] = n + 1
    end
    expected[1, end] = n + 1
    @test evaluate.(result.matrix) == expected
    @test result.grading == Diagonal([(2 * i - n)//2 for i in 0:n])
    @test result.metric^2 == Matrix{QQ}(I, n + 1, n + 1)
    @test fano_curve_classes(X) == (n + 1, [[0], [1]])
  end
  P1 = projective_space(1)
  @test small_quantum_multiplication(P1; betas=[[0], [1], [1]]).matrix ==
    small_quantum_multiplication(P1).matrix
  @test evaluate.(small_quantum_multiplication(P1; betas=[[0]]).matrix) == [0 0; 2 0]

  F = full_flag_variety(TypeA{2})
  @test evaluate.(small_quantum_multiplication(F).matrix) ==
    [0 2 2 0 0 4; 2 0 0 2 0 0; 2 0 0 0 2 0; 0 2 4 0 0 2; 0 4 2 0 0 2; 0 0 0 2 2 0]

  X = projective_space(4)
  quadric = small_quantum_multiplication(X; bundle=line_bundle(X, 2))
  @test evaluate.(quadric.matrix) == [0 0 6 0; 3 0 0 6; 0 3 0 0; 0 0 3 0]
  G = Gr(2, 4)
  quotient = small_quantum_multiplication(G; bundle=universal_quotient_bundle(G))
  @test evaluate.(quotient.matrix) == [0 0 3; 3 0 0; 0 3 0]

  S = product(P1, P1)
  E = line_bundle(S, [1, 2])
  @test_throws ArgumentError fano_curve_classes(S; bundle=E)
  @test_throws ArgumentError small_quantum_multiplication(S; bundle=E)
  result = small_quantum_multiplication(zero_locus(E); betas=[[0, 0], [2, 1]])
  @test evaluate.(result.matrix) == [0 4; 1 0]
  K3 = zero_locus(line_bundle(projective_space(3), 4))
  @test fano_curve_classes(K3) == (0, [[0]])
  @test all(iszero, small_quantum_multiplication(K3).matrix)
  Y = projective_space(2)
  @test_throws ArgumentError fano_curve_classes(Y; bundle=line_bundle(Y, 4))
  @test size(
    small_quantum_multiplication(Y; bundle=line_bundle(Y, 4), betas=[[0]]).matrix
  ) == (2, 2)
end

@testset "Nonlinear defining bundles through degree three" begin
  # These low-degree parts of the full fourfold fixtures exercise exterior
  # and symmetric powers without the much larger degree-five tree sums.
  G = Gr(3, 7)
  E =
    exterior_power(dual(universal_subbundle(G)), 2) +
    exterior_power(universal_quotient_bundle(G), 3) + line_bundle(G, 1)
  c5 = small_quantum_multiplication(zero_locus(E); betas=[[0], [1], [2], [3]])
  q1 = NovikovPolynomial{1}(Dict((1,) => 1))
  z, u = zero(q1), one(q1)
  q2 = q1 * q1
  q3 = q2 * q1
  @test c5.matrix == [
    z 72 * q2 732 * q3 888 * q3 z z
    u 5 * q1 92 * q2 112 * q2 324 * q3 z
    z u z 4 * q1 12 * q2 108 * q3
    z u 6 * q1 4 * q1 24 * q2 180 * q3
    z z 5 * u 6 * u 5 * q1 72 * q2
    z z z z u z
  ]

  G = Gr(2, 7)
  E = 3 * line_bundle(G, 1) + symmetric_power(dual(universal_subbundle(G)), 2)
  b8 = small_quantum_multiplication(zero_locus(E); betas=[[0], [1], [2], [3]])
  @test b8.matrix == [
    z 96 * q2 936 * q3 1656 * q3 z z
    u 6 * q1 108 * q2 192 * q2 2160 * q3 z
    z u 2 * q1 12 * q1 120 * q2 648 * q3
    z u 4 * q1 2 * q1 100 * q2 504 * q3
    z z u (9//5) * u 6 * q1 (288//5) * q2
    z z z z (5//3) * u z
  ]
end
