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
  @test string.(result.matrix) == [
    "0" "72//1*q1^2" "732//1*q1^3" "888//1*q1^3" "2496//1*q1^4" "18720//1*q1^5"
    "1//1" "5//1*q1" "92//1*q1^2" "112//1*q1^2" "324//1*q1^3" "2496//1*q1^4"
    "0" "1//1" "0" "4//1*q1" "12//1*q1^2" "108//1*q1^3"
    "0" "1//1" "6//1*q1" "4//1*q1" "24//1*q1^2" "180//1*q1^3"
    "0" "0" "5//1" "6//1" "5//1*q1" "72//1*q1^2"
    "0" "0" "0" "0" "1//1" "0"
  ]
  @info "Fourfold regression completed: Küchle c5"
  flush(stderr)
end
