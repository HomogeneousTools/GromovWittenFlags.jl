@testset "Recursive branch localization" begin
  P1, P2, P3 = projective_space(1), projective_space(2), projective_space(3)
  product_space = product(P1, P1)
  flag = full_flag_variety(TypeA{2})
  grass = Gr(2, 4)
  cases = (
    (P1, [], [1], nothing, 1),
    (P1, [point_class(P1)], [1], nothing, 1),
    (P1, fill(point_class(P1), 2), [1], nothing, 1),
    (P2, fill(point_class(P2), 2), [1], nothing, 1),
    (P2, fill(point_class(P2), 5), [2], nothing, 1),
    (P2, fill(point_class(P2), 8), [3], nothing, 12),
    (product_space, fill(point_class(product_space), 3), [1, 1], nothing, 1),
    (product_space, fill(point_class(product_space), 5), [2, 1], nothing, 1),
    (flag, fill(point_class(flag), 2), [1, 1], nothing, 1),
    (grass,
      [schubert_class(grass, [2]), schubert_class(grass, [1, 3, 2]),
        point_class(grass)], [1], nothing, 1),
    (P3, [], [1], line_bundle(P3, 3), 27),
  )
  for (X, insertions, beta, bundle, answer) in cases
    @test gromov_witten_invariant(X, insertions, beta; bundle, algorithm=:branches) ==
      answer
    @test gromov_witten_invariant(X, insertions, beta; bundle, algorithm=:branches) ==
      gromov_witten_invariant(X, insertions, beta; bundle)
  end

  # Both arithmetic implementations use the same recursion, compared with
  # the independent enumerator at individual torus points, not only after summation.
  for X in (P2, full_flag_variety(TypeB{2}), full_flag_variety(TypeG2))
    data = GWF._gw_data(X)
    beta = X == P2 ? [2] : [1, 0]
    classes = [1, length(data.representatives)]
    for T in (QQ, GWF.QQFieldElem)
      evaluation = GWF._gw_evaluation(data, 0, T)
      geometry_before = deepcopy((evaluation.vertex_euler, evaluation.root_values))
      expected = only(GWF._positive_degree_invariants((evaluation,), classes, beta))
      @test GWF._branch_invariant(evaluation, classes, beta) == expected
      @test GWF._branch_invariant(evaluation, classes, beta) == expected
      @test (evaluation.vertex_euler, evaluation.root_values) == geometry_before
    end
  end

  # The nonsplit presentation of P2 forces H0-H1 cancellation on double covers.
  F = flag_variety(4, [1, 3])
  L, B = tautological_bundles(F)
  twist = ConvexEulerTwist(dual(L) + dual(B);
    certificate="U3 dual is globally generated, and its zero locus is P2")
  reference = small_quantum_multiplication(F; bundle=twist, betas=[[0, 0], [1, 0]])
  recursive = small_quantum_multiplication(F; bundle=twist,
    betas=[[0, 0], [1, 0]], algorithm=:branches)
  @test recursive.matrix == reference.matrix
  @test recursive.metric == reference.metric
  @test gromov_witten_invariant(F, fill(reference.basis[3], 5), [2, 0];
    bundle=twist, algorithm=:branches) == 1

  for (X, bundle) in ((P1, nothing), (P2, nothing), (grass, nothing),
    (P3, line_bundle(P3, 2)), (P3, line_bundle(P3, 3)))
    tree = small_quantum_multiplication(X; bundle)
    branch = small_quantum_multiplication(X; bundle, algorithm=:branches)
    @test branch.matrix == tree.matrix
    @test branch.basis == tree.basis
    @test branch.metric == tree.metric
  end

  @test gromov_witten_invariant(P1, fill(point_class(P1), 2), [0]; algorithm=:branches) == 0
  @test gromov_witten_invariant(P1, [point_class(P1), first(schubert_classes(P1))],
    [0]; algorithm=:branches) == 1
  @test_throws ArgumentError gromov_witten_invariant(P1, [], [1]; algorithm=:unknown)
  @test_throws ArgumentError small_quantum_multiplication(P1; algorithm=:unknown)
  @test_throws ArgumentError GWF._branch_mask_splits(8)
  events = Any[]
  @test gromov_witten_invariant(P2, fill(point_class(P2), 2), [1];
    algorithm=:branches, progress=event -> push!(events, event)) == 1
  @test any(event -> event.phase == :geometry, events)
  @test count(event -> event.phase == :complete, events) == 2
end
