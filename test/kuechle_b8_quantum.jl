@testset "Küchle b8 quantum multiplication" begin
  @info "Fourfold regression: Küchle b8"
  flush(stderr)
  G = Gr(2, 7)
  S = universal_subbundle(G)
  bundle = 3 * line_bundle(G, 1) + symmetric_power(dual(S), 2)
  X = zero_locus(bundle)

  @test dimension(X) == 4

  result = small_quantum_multiplication(X)
  @test string.(result.matrix) == [
    "0" "96//1*q1^2" "936//1*q1^3" "1656//1*q1^3" "18560//1*q1^4" "93120//1*q1^5"
    "1//1" "6//1*q1" "108//1*q1^2" "192//1*q1^2" "2160//1*q1^3" "11136//1*q1^4"
    "0" "1//1" "2//1*q1" "12//1*q1" "120//1*q1^2" "648//1*q1^3"
    "0" "1//1" "4//1*q1" "2//1*q1" "100//1*q1^2" "504//1*q1^3"
    "0" "0" "1//1" "9//5" "6//1*q1" "288//5*q1^2"
    "0" "0" "0" "0" "5//3" "0"
  ]
  @info "Fourfold regression completed: Küchle b8"
  flush(stderr)
end
