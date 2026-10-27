@testset "Exact Novikov polynomials" begin
  P = NovikovPolynomial{2}
  p = P(Dict((1, 0) => 2, (0, 2) => -1//3, (0, 0) => 0))
  q = P(Dict((0, 1) => 3//2))
  @test iszero(zero(P))
  @test iszero(p - p)
  @test iszero(0 * p)
  @test one(p) * p == p
  @test zero(p) + p == p
  @test p + q == q + p
  @test p * q == q * p
  @test evaluate(p, [2, 3]) == 1
  @test evaluate(p + q, [2, 3]) == 11//2
  @test evaluate(p * q, [2, 3]) == 9//2
  @test evaluate(3//2 * p, [2, 3]) == 3//2
  @test evaluate(p * 2) == 10//3
  @test evaluate(p) isa QQ
  @test P(2) == P(Dict((0, 0) => 2))
  @test hash(copy(p)) == hash(p)
  @test isequal(copy(p), p)
  @test copy(p).terms !== p.terms
  @test string(p) == "-1//3*q2^2 + 2//1*q1"
  @test string(-q) == "-3//2*q2"
  @test string(zero(P)) == "0"
  @test evaluate(NovikovPolynomial{0}(2), []) == 2
  @test_throws ArgumentError P(Dict((-1, 0) => 1))
  @test_throws ArgumentError P(Dict((1,) => 1))
  @test_throws ArgumentError P(Dict((1, 0) => 1.5))
  @test_throws ArgumentError P(Dict((1, 0) => 1//0))
  @test_throws ArgumentError P(1//0)
  @test_throws ArgumentError (1//0) * p
  @test_throws ArgumentError (1//0) * zero(P)
  @test_throws ArgumentError P(Dict((big(typemax(Int)) + 1, 0) => 1))
  @test_throws ArgumentError NovikovPolynomial{-1}()
  @test_throws ArgumentError NovikovPolynomial{-1}(2)
  @test_throws ArgumentError NovikovPolynomial{2.0}()
  @test_throws ArgumentError evaluate(p, [1])
  @test_throws OverflowError P(Dict((typemax(Int), 0) => 1)) * p
  @test P(Dict{Any,Int}((1, 0) => 2, [1, 0] => 3)) == P(Dict((1, 0) => 5))
  @test iszero(P(Dict{Any,Int}((1, 0) => 2, [1, 0] => -2)))
  @test evaluate.(p) == evaluate(p)
  @test evaluate(p, (2, 3)) == 1
  @test evaluate(P(Dict((1, 0) => 1)), [0.5, 1]) == 1//2
  @test_throws ArgumentError evaluate(p, [Inf, 1])
  @test_throws ArgumentError evaluate(p, [NaN, 1])
  @test_throws ArgumentError evaluate(p, [1//0, 1])
  @test_throws ArgumentError evaluate(p, "11")
  @test zero(P) != zero(NovikovPolynomial{1})
  @test !isequal(zero(P), zero(NovikovPolynomial{1}))
  @test_throws DimensionMismatch p + NovikovPolynomial{1}(1)
  @test_throws DimensionMismatch p - NovikovPolynomial{1}(1)
  @test_throws DimensionMismatch p * NovikovPolynomial{1}(1)
  source = Dict((1, 0) => 1)
  owned = P(source)
  empty!(source)
  @test owned == P(Dict((1, 0) => 1))
end

@testset "Streaming balanced exact sums" begin
  accumulator = GWF._ExactSum(QQ)
  @test GWF._sum_value(accumulator) == 0
  expected = zero(QQ)
  for i in 1:1025
    value = big(mod(i, 37) - 18) // big(mod(i, 29) + 1)
    GWF._add_exact!(accumulator, value)
    expected += value
    @test GWF._sum_value(accumulator) == expected
  end
  @test length(accumulator.partials) <= 11
  count = accumulator.count
  GWF._add_exact!(accumulator, zero(QQ))
  @test accumulator.count == count
  GWF._add_exact!(accumulator, -expected)
  @test iszero(GWF._sum_value(accumulator))
end

@testset "Exact rank and radical quotient" begin
  A = QQ[1 2 3; 2 4 6; 0 1 1]
  @test GWF._exact_rank(A) == 2
  @test GWF._independent_indices(A) == [1, 3]
  @test GWF._exact_rank(zeros(QQ, 0, 0)) == 0
  @test GWF._exact_rank(zeros(QQ, 3, 2)) == 0
  @test A == QQ[1 2 3; 2 4 6; 0 1 1]
  @test GWF._independent_indices(QQ[0 0 1; 1 0 2; 0 1 3; 1 1 6]) == [1, 2, 3]
  @test GWF._independent_indices(QQ[0 0 1; 0 0 2; 1 2 0]) == [1, 3]
  @test GWF._exact_rank(transpose(A)) == 2
  @test GWF._exact_rank(@view A[1:2, :]) == 1
  @test GWF._independent_indices(Rational{Int}[1 2; 2 4; 1 3]) == [1, 3]
  @test GWF._exact_rank(zeros(QQ, 3, 0)) == 0
end
