@testset "FLINT quantum blocks preserve exact rational coefficients" begin
  P3 = projective_space(3)
  examples = ((projective_space(2), [1], nothing),
    (Gr(2, 4), [1], nothing),
    (full_flag_variety(TypeA{2}), [1, 1], nothing),
    (P3, [1], line_bundle(P3, 3)))
  for (X, beta, bundle) in examples
    data = GWF._gw_data(X)
    dim, c_full, weights = GWF._complete_intersection_data(data, bundle)
    c = [c_full[node] for node in marked_nodes(X)]
    expected_degree = dim - 1 + sum(c .* beta)
    classes = collect(eachindex(data.representatives))
    # Only top-degree insertions give scalar, torus-independent integrals.
    pairs = [
      (i, j) for i in classes for j in i:length(classes)
      if data.codimensions[i] + data.codimensions[j] == expected_degree
    ]
    @test !isempty(pairs)
    reference = GWF._checked_torus_evaluations(data, weights) do evaluations
      foreach(evaluation -> GWF._prepare_edges!(evaluation, beta), evaluations)
      GWF._positive_degree_pair_blocks(evaluations, classes, pairs, beta)
    end
    actual = GWF._checked_pair_block(data, classes, pairs, beta, weights)
    @test eltype(actual) == QQ
    @test actual == reference
    @test actual == GWF._checked_pair_block(data, classes, pairs, beta, weights)
  end
end

@testset "FLINT evaluations retain the pole and disagreement guards" begin
  T = GWF.QQFieldElem
  data = GWF._gw_data(projective_space(1))
  @test_throws GWF._LocalizationPole GWF._divide_or_pole(one(T), zero(T))
  calls = Ref(0)
  value = GWF._checked_torus_evaluations(data, nothing, T) do evaluations
    @test eltype(first(evaluations).point) == T
    calls[] += 1
    calls[] == 1 && throw(GWF._LocalizationPole())
    (T(7), T(7))
  end
  @test calls[] == 2
  @test value == 7
  @test_throws ErrorException GWF._checked_torus_evaluations(data, nothing, T) do _
    (one(T), T(2))
  end
  @test_throws DivideError GWF._checked_torus_evaluations(data, nothing, T) do _
    div(1, 0)
  end
end
