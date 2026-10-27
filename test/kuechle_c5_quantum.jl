@testset "Example 6.3: Kuchle c5 fourfold" begin
  @info "Fourfold regression: Küchle c5"
  flush(stderr)
  X = Gr(3, 7)
  bundle = reduce(
    direct_sum,
    [
      exterior_power(dual(universal_subbundle(X)), 2),
      exterior_power(universal_quotient_bundle(X), 3),
      line_bundle(X, 1),
    ],
  )
  result = small_quantum_multiplication(zero_locus(bundle))
  q1 = NovikovPolynomial{1}(Dict((1,) => 1))
  z, u = zero(q1), one(q1)
  q2 = q1 * q1
  q3 = q2 * q1
  q4 = q3 * q1
  q5 = q4 * q1
  @test result.matrix == [
    z 72 * q2 732 * q3 888 * q3 2496 * q4 18720 * q5
    u 5 * q1 92 * q2 112 * q2 324 * q3 2496 * q4
    z u z 4 * q1 12 * q2 108 * q3
    z u 6 * q1 4 * q1 24 * q2 180 * q3
    z z 5 * u 6 * u 5 * q1 72 * q2
    z z z z u z
  ]
  @info "Fourfold regression completed: Küchle c5"
  flush(stderr)
end
