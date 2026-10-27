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
  @test string.(c5.matrix) == [
    "0" "72//1*q1^2" "732//1*q1^3" "888//1*q1^3" "0" "0"
    "1//1" "5//1*q1" "92//1*q1^2" "112//1*q1^2" "324//1*q1^3" "0"
    "0" "1//1" "0" "4//1*q1" "12//1*q1^2" "108//1*q1^3"
    "0" "1//1" "6//1*q1" "4//1*q1" "24//1*q1^2" "180//1*q1^3"
    "0" "0" "5//1" "6//1" "5//1*q1" "72//1*q1^2"
    "0" "0" "0" "0" "1//1" "0"
  ]

  G = Gr(2, 7)
  E = 3 * line_bundle(G, 1) + symmetric_power(dual(universal_subbundle(G)), 2)
  b8 = small_quantum_multiplication(zero_locus(E); betas=[[0], [1], [2], [3]])
  @test string.(b8.matrix) == [
    "0" "96//1*q1^2" "936//1*q1^3" "1656//1*q1^3" "0" "0"
    "1//1" "6//1*q1" "108//1*q1^2" "192//1*q1^2" "2160//1*q1^3" "0"
    "0" "1//1" "2//1*q1" "12//1*q1" "120//1*q1^2" "648//1*q1^3"
    "0" "1//1" "4//1*q1" "2//1*q1" "100//1*q1^2" "504//1*q1^3"
    "0" "0" "1//1" "9//5" "6//1*q1" "288//5*q1^2"
    "0" "0" "0" "0" "5//3" "0"
  ]
end
