@testset "Exact scalar local monodromy" begin
  import Nemo
  R, t = Nemo.polynomial_ring(Nemo.QQ, :t)
  sizes(p, point=:zero) = [b.sizes for b in local_monodromy(p, point).blocks]

  # Resonant exponents 0,1 allow both apparent and logarithmic behavior.
  @test sizes([[0], [-1], [1]]) == [[1, 1]]
  @test sizes([[0, -1], [-1], [1]]) == [[2]]
  @test sizes([[0], [0], [1]]) == [[2]]
  @test sizes([[0], [0], [0], [0], [1]]) == [[4]]
  @test sizes([[0], [2], [1]]) == [[1, 1]] # negative exponent -2
  @test_throws ArgumentError local_monodromy([[0], [-1000], [1]])

  # The nonunipotent part contributes to rank(T-I).
  fractional = ramification_data([[-1//2], [1]])
  @test fractional.ramification == 2
  @test fractional.defect == 0
  @test all(p.fixed_dimension == 0 for p in fractional.local_data)
  @test log_monodromy_matrix(fractional.local_data[1]) == Nemo.matrix(Nemo.QQ, 1, 1, [1//2])
  @test ramification_defect(fractional) == 0
  @test ramification_defect([[-1//2], [1]]) == 0
  # A trivial/reducible local system need not obey the irreducible lower bound.
  @test ramification_data([[0], [-1], [1]]).defect == -4

  # A simple zero with the same indicial roots 0,2 can be apparent or logarithmic.
  @test sizes([[0], [1, -2], [-1, 1]], t-1) == [[1, 1]]
  @test sizes([[0, 0, 1], [1, -2], [-1, 1]], t-1) == [[2]]
  # A repeated leading zero exercises the general number-field Frobenius path.
  @test sizes([[0], [-1, 1], [1, -2, 1]], t-1) == [[2]]
  f = t^2-2
  @test sizes([4t^4, -2t^2*f-2f^2, f^2], f) == [[2]]
  ordinary = local_monodromy([[0], [0], [1]], t-1)
  @test ordinary.apparent
  @test ordinary.ramification == 0
  @test local_monodromy([[0, 0], [0, 1]], :zero).apparent

  # Appendix B.1 of arXiv:1406.4891: P4 has five conifolds and MUM at zero.
  p4 = [[0, 0, 0, 0, 0, 75000], [0, 0, 0, 0, 0, 156250],
    [0, 0, 0, 0, 0, 109375], [0, 0, 0, 0, 0, 31250],
    [-1, 0, 0, 0, 0, 3125]]
  data = ramification_data(p4)
  @test data.order == 4
  @test data.ramification == 8
  @test data.defect_zero
  @test data.local_data[1].blocks[1].sizes == [4]
  @test sum(p.degree for p in data.local_data[2:(end - 1)]) == 5
  @test all(p.ramification == 1 for p in data.local_data[2:(end - 1)])
  @test data.local_data[end].apparent
  @test Nemo.rank(log_monodromy_matrix(data.local_data[1])) == 3

  # Appendix B.8: P2 x P2 has six conifolds, giving defect one.
  product = [[0, 0, 0, 972, 0, 0, 233280], [0, 0, 0, 2592, 0, 0, 454896],
    [0, 0, 0, 2565, 0, 0, 285768], [0, 0, 0, 1134, 0, 0, 69984],
    [-1, 0, 0, 189, 0, 0, 5832]]
  @test ramification_data(product).ramification == 9
  @test ramification_defect(product) == 1

  # Coefficient content and a common polynomial factor do not change the local system.
  polynomials = [R(row) for row in p4]
  @test ramification_defect([(t-2)*p for p in polynomials]) == 0
  @test_throws ArgumentError local_monodromy([[1.0], [1]])
  @test_throws ArgumentError local_monodromy([[1]])
  @test_throws ArgumentError local_monodromy([[0], [0]])
  @test_throws ArgumentError local_monodromy([[1], [0], [0, 1]])
  @test_throws ArgumentError local_monodromy([[-2], [0], [1]])
  @test_throws ArgumentError local_monodromy(p4, (t-1)*(t-2))
  @test_throws ArgumentError local_monodromy(p4, :unknown)
end

@testset "arXiv:1406.4891 Appendix B (35 operators)" begin
  import Nemo
  import TOML
  fixture = TOML.parsefile(joinpath(@__DIR__, "data", "monodromy-1406.4891.toml"))
  R, _ = Nemo.polynomial_ring(Nemo.QQ, :t)
  rational(s) = begin
    parts = split(s, '/')
    Nemo.QQ(
      parse(BigInt, parts[1]), length(parts) == 1 ? BigInt(1) : parse(BigInt, parts[2])
    )
  end
  @test length(fixture["cases"]) == 35
  @test count(c -> c["expected_defect"] == 0, fixture["cases"]) == 24
  for case in fixture["cases"]
    @testset "$(case["appendix"])" begin
      coefficients = [rational.(row) for row in case["coefficients"]]
      data = ramification_data(coefficients)
      @test data.defect == case["expected_defect"]
      for point in case["expected_local"]
        row = point["matrix"]
        n = length(row)
        expected = Nemo.matrix(Nemo.QQ, [rational(row[i][j]) for i in 1:n, j in 1:n])
        loci = if point["locus"] == "finite"
          [f for (f, _) in Nemo.factor(R(rational.(point["locus_coefficients"])))]
        else
          [Symbol(point["locus"])]
        end
        for locus in loci
          actual = log_monodromy_matrix(local_monodromy(coefficients, locus))
          @test size(actual) == size(expected)
          for exponent in unique([expected[i, i] for i in 1:n]), power in 1:n
            @test Nemo.rank((actual - exponent * one(actual))^power) ==
              Nemo.rank((expected - exponent * one(expected))^power)
          end
        end
      end
    end
  end
end
