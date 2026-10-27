@testset "Küchle b8 quantum multiplication" begin
  @info "Fourfold regression: Küchle b8"
  flush(stderr)
  G = Gr(2, 7)
  S = universal_subbundle(G)
  bundle = 3 * line_bundle(G, 1) + symmetric_power(dual(S), 2)
  X = zero_locus(bundle)

  @test dimension(X) == 4

  result = small_quantum_multiplication(X)
  q1 = NovikovPolynomial{1}(Dict((1,) => 1))
  z, u = zero(q1), one(q1)
  q2 = q1 * q1
  q3 = q2 * q1
  q4 = q3 * q1
  q5 = q4 * q1
  @test result.matrix == [
    z 96 * q2 936 * q3 1656 * q3 18560 * q4 93120 * q5
    u 6 * q1 108 * q2 192 * q2 2160 * q3 11136 * q4
    z u 2 * q1 12 * q1 120 * q2 648 * q3
    z u 4 * q1 2 * q1 100 * q2 504 * q3
    z z u (9//5) * u 6 * q1 (288//5) * q2
    z z z z (5//3) * u z
  ]
  @info "Fourfold regression completed: Küchle b8"
  flush(stderr)
end
