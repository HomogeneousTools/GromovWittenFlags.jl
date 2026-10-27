@testset "arXiv:2609.21914 examples" begin
  matrix_strings(result) = string.(result.matrix)

  @testset "Example 6.1: cubic fourfold" begin
    @info "Fourfold regression: cubic"
    X = projective_space(5)
    result = small_quantum_multiplication(zero_locus(line_bundle(X, 3)))
    @test matrix_strings(result) == [
      "0" "0" "18//1*q1" "0" "0"
      "3//1" "0" "0" "45//1*q1" "0"
      "0" "3//1" "0" "0" "18//1*q1"
      "0" "0" "3//1" "0" "0"
      "0" "0" "0" "3//1" "0"
    ]
  end

  @testset "Example 6.5: ordinary Gushel--Mukai fourfold" begin
    @info "Fourfold regression: Gushel–Mukai"
    X = Gr(2, 5)
    bundle = line_bundle(X, 1) + line_bundle(X, 2)
    result = small_quantum_multiplication(zero_locus(bundle))
    @test matrix_strings(result) == [
      "0" "12//1*q1" "0" "0" "96//1*q1^2" "0"
      "2//1" "0" "12//1*q1" "20//1*q1" "0" "48//1*q1^2"
      "0" "2//1" "0" "0" "8//1*q1" "0"
      "0" "2//1" "0" "0" "16//1*q1" "0"
      "0" "0" "2//1" "3//1" "0" "6//1*q1"
      "0" "0" "0" "0" "4//1" "0"
    ]
  end

  @testset "Example 6.7: quartic fourfold" begin
    @info "Fourfold regression: quartic"
    X = projective_space(5)
    result = small_quantum_multiplication(zero_locus(line_bundle(X, 4)))
    @test matrix_strings(result) == [
      "0" "48//1*q1" "0" "5568//1*q1^2" "0"
      "2//1" "0" "208//1*q1" "0" "5568//1*q1^2"
      "0" "2//1" "0" "208//1*q1" "0"
      "0" "0" "2//1" "0" "48//1*q1"
      "0" "0" "0" "2//1" "0"
    ]
  end

  @testset "Example 6.8: quadric--cubic fourfold" begin
    @info "Fourfold regression: quadric–cubic"
    X = projective_space(6)
    bundle = line_bundle(X, 2) + line_bundle(X, 3)
    result = small_quantum_multiplication(zero_locus(bundle))
    @test matrix_strings(result) == [
      "0" "24//1*q1" "0" "1152//1*q1^2" "0"
      "2//1" "0" "84//1*q1" "0" "1152//1*q1^2"
      "0" "2//1" "0" "84//1*q1" "0"
      "0" "0" "2//1" "0" "24//1*q1"
      "0" "0" "0" "2//1" "0"
    ]
  end

  @testset "Example 6.9: intersection of three quadrics" begin
    @info "Fourfold regression: three quadrics"
    X = projective_space(7)
    result = small_quantum_multiplication(zero_locus(3 * line_bundle(X, 2)))
    @test matrix_strings(result) == [
      "0" "16//1*q1" "0" "448//1*q1^2" "0"
      "2//1" "0" "48//1*q1" "0" "448//1*q1^2"
      "0" "2//1" "0" "48//1*q1" "0"
      "0" "0" "2//1" "0" "16//1*q1"
      "0" "0" "0" "2//1" "0"
    ]
  end

  @testset "Example 6.11: K3-29" begin
    @info "Fourfold regression: K3-29"
    X = product(projective_space(2), projective_space(4))
    bundle = line_bundle(X, [1, 2]) + line_bundle(X, [1, 1])
    result = small_quantum_multiplication(zero_locus(bundle))
    @test matrix_strings(result) == [
      "0" "0" "4//1*q2" "30//1*q1*q2" "12//1*q1*q2" "30//1*q1*q2" "56//1*q1^2*q2" "32//1*q1^2*q2" "60//1*q1*q2^2 + 15//1*q1^3*q2"
      "1//1" "-q1" "0" "4//1*q2" "0" "10//1*q2" "33//1*q1*q2" "18//1*q1*q2" "12//1*q1^2*q2"
      "2//1" "3//1*q1" "0" "0" "0" "4//1*q2" "54//1*q1*q2" "24//1*q1*q2" "48//1*q1^2*q2"
      "0" "2//1" "1//1" "-q1" "0" "0" "4//1*q2" "0" "18//1*q1*q2"
      "0" "1//1" "0" "0" "-q1" "0" "10//1*q2" "4//1*q2" "9//1*q1*q2"
      "0" "0" "2//1" "3//1*q1" "2//1*q1" "0" "0" "0" "18//1*q1*q2"
      "0" "0" "0" "2//1" "0" "4//1" "7//2*q1" "3//1*q1" "0"
      "0" "0" "0" "1//1" "2//1" "-7//2" "-21//4*q1" "-9//2*q1" "6//1*q2"
      "0" "0" "0" "0" "0" "0" "8//3" "4//3" "0"
    ]
  end
end
