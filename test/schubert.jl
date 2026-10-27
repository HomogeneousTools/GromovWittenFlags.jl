@testset "Schubert basis and curve coordinates" begin
  X = projective_space(2)
  basis = schubert_classes(X)
  @test length.(basis) == 0:2
  @test point_class(X) == last(basis)
  @test schubert_class(X, ()) == first(basis)
  @test schubert_class(X, (1,)) == schubert_class(X, [1])
  pop!(basis)
  @test length(schubert_classes(X)) == 3

  G = Gr(2, 4)
  data = GWF._gw_data(G)
  @test length.(schubert_classes(G)) == [0, 1, 2, 2, 3, 4]
  @test GWF._normalize_beta(data, [1]) == [1]
  @test GWF._normalize_beta(data, (0, 1, 0)) == [1]
  @test_throws ArgumentError GWF._normalize_beta(data, [1, 0, 0])
  @test_throws ArgumentError GWF._normalize_beta(data, [1, 0])
  @test_throws ArgumentError GWF._normalize_beta(data, [-1])
  @test_throws ArgumentError GWF._normalize_beta(data, [1.0])
  @test_throws ArgumentError GWF._normalize_beta(data, "1")
  @test_throws ArgumentError GWF._normalize_beta(data, [big(typemax(Int)) + 1])
  @test_throws ArgumentError schubert_class(G, [0])
  @test_throws ArgumentError schubert_class(G, [4])

  # Both the collection and each reduced word must be owned at the boundary.
  classes = schubert_classes(G)
  empty!(GWF.word(last(classes)))
  empty!(GWF.word(point_class(G)))
  empty!(GWF.word(schubert_class(G, [2])))
  @test length(point_class(G)) == 4
  @test length(schubert_class(G, [2])) == 1
  @test length.(schubert_classes(G)) == [0, 1, 2, 2, 3, 4]

  for X in (full_flag_variety(TypeB{2}), full_flag_variety(TypeG2),
    product(projective_space(1), projective_space(2)))
    @test length(schubert_classes(X)) == euler_characteristic(X)
    @test length(point_class(X)) == dimension(X)
  end
  @test clear_localization_caches!() === nothing
  @test isempty(GWF._GW_DATA_CACHE)
  @test length(schubert_classes(G)) == 6
end
