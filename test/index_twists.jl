@testset "Certified convex index twists" begin
  P3 = projective_space(3)
  E = line_bundle(P3, 3)
  @test_throws ArgumentError ConvexEulerTwist(E; certificate=" ")
  ordinary = small_quantum_multiplication(P3; bundle=E, betas=[[0], [1]])
  certified = small_quantum_multiplication(P3;
    bundle=ConvexEulerTwist(E; certificate="O(3) is globally generated"),
    betas=[[0], [1]])
  @test ordinary.matrix == certified.matrix
  @test ordinary.metric == certified.metric

  # The standard SL2 representation is a trivial rank-two bundle with
  # graded class O(1)+O(-1). On even covers its intermediate H0 and H1
  # characters are zero and must cancel before evaluating an Euler factor.
  P1 = projective_space(1)
  trivial_graded = line_bundle(P1, 1) + line_bundle(P1, -1)
  data1 = GWF._gw_data(P1)
  indexed = GWF._gw_evaluation(data1, 0; weights=GWF._bundle_weights(trivial_graded))
  untwisted = GWF._gw_evaluation(data1, 0)
  for cover in 1:5, color in 1:2
    @test GWF._edge_factor(indexed, color, 1, cover) ==
      GWF._edge_factor(untwisted, color, 1, cover) *
          prod(indexed.fiber_values[color, :])
  end

  # U3 dual on Fl(1,3;4) is globally generated; its regular zero locus is P2.
  # Its graded class has negative splittings and H1 on double covers.
  F = flag_variety(4, [1, 3])
  L, B = tautological_bundles(F)
  graded = dual(L) + dual(B)
  twist = ConvexEulerTwist(graded;
    certificate="U3 dual is globally generated; a regular section cuts out P2")
  @test_throws ArgumentError small_quantum_multiplication(F;
    bundle=graded, betas=[[0, 0]])
  result = small_quantum_multiplication(F; bundle=twist, betas=[[0, 0], [1, 0]])
  @test result.metric == [0 0 1; 0 1 0; 1 0 0]
  @test evaluate.(result.matrix) == [0 0 3; 3 0 0; 0 3 0]
  @test evaluate.(result.matrix, Ref((2, 7))) == [0 0 6; 3 0 0; 0 3 0]
  @test gromov_witten_invariant(F, fill(result.basis[3], 2), [1, 0]; bundle=twist) == 1
  @test gromov_witten_invariant(F, fill(result.basis[3], 5), [2, 0]; bundle=twist) == 1

  # Independently, U3 is constant along a (1,0) curve, so its dual
  # restricts to O^3 and its index has rank three on every cover.
  data = GWF._gw_data(F)
  _, _, weights = GWF._complete_intersection_data(data, twist)
  evaluation = GWF._gw_evaluation(data, 0; weights)
  roots = findall(==([1, 0]), data.curve_classes)
  @test !isempty(roots)
  for root in roots, cover in 1:5
    terms = GWF._build_edge_index_terms(data, evaluation.weight_coordinates,
      evaluation.splittings, root, cover)
    @test sum(last, terms) == 3
    @test all(term -> last(term) > 0, terms)
    for color in eachindex(data.representatives)
      endpoint = data.steps[color, root]
      reverse = data.edge_roots[endpoint][color]
      @test GWF._edge_factor(evaluation, color, root, cover) ==
        GWF._edge_factor(evaluation, endpoint, reverse, cover)
    end
  end

  # A false convexity assertion cannot hide a negative character index.
  P2 = projective_space(2)
  bad = ConvexEulerTwist(line_bundle(P2, -1);
    certificate="deliberately false certificate for the rejection test")
  @test_throws ArgumentError gromov_witten_invariant(P2,
    fill(point_class(P2), 6), [2]; bundle=bad)
end
