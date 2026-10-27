@testset "Genus-zero primary invariants" begin
  X = projective_space(2)
  pt = point_class(X)
  unit, h = schubert_classes(X)[1:2]
  @test expected_gromov_witten_degree(X, [1], 2) == 4
  @test expected_gromov_witten_degree(X, [0], 2) == 2
  @test expected_gromov_witten_degree(X, [0], 4) == 3
  @test_throws ArgumentError expected_gromov_witten_degree(X, [1], -1)
  @test gromov_witten_invariant(X, [pt, pt], [1]) == 1
  @test gromov_witten_invariant(X, [pt, pt], [1]) isa QQ
  @test gromov_witten_invariant(X, [pt, pt, h], [1]) == 1
  @test gromov_witten_invariant(X, fill(pt, 5), [2]) == 1
  @test gromov_witten_invariant(X, [pt, unit], [0]) == 1
  @test gromov_witten_invariant(X, [h, h, unit], [0]) == 1
  @test gromov_witten_invariant(X, [pt, unit, unit, unit], [0]) == 0
  @test gromov_witten_invariant(X, [pt, h, unit, unit], [0]) == 0
  @test gromov_witten_invariant(X, [h, h], [1]) == 0
  @test_throws ArgumentError gromov_witten_invariant(X, [pt, pt], [-1])
  @test_throws ArgumentError gromov_witten_invariant(
    X, [point_class(projective_space(3))], [1]
  )
  @test_throws ArgumentError gromov_witten_invariant(X, [1], [1])

  G = Gr(2, 4)
  s1 = schubert_class(G, [2])
  s21 = schubert_class(G, [1, 3, 2])
  classes = [s1, s21, point_class(G)]
  @test gromov_witten_invariant(G, classes, [1]) == 1
  @test gromov_witten_invariant(G, classes, [0, 1, 0]) == 1
  F = full_flag_variety(TypeA{2})
  @test gromov_witten_invariant(F, fill(point_class(F), 2), [1, 1]) == 1
  B = full_flag_variety(TypeB{2})
  @test gromov_witten_invariant(B, [schubert_class(B, [1]), point_class(B)], [1, 0]) == 1

  P3 = projective_space(3)
  E = line_bundle(P3, 3)
  Z = zero_locus(E)
  @test gromov_witten_invariant(P3, [], [1]; bundle=E) == 27
  @test gromov_witten_invariant(Z, [], [1]) == 27
  @test expected_gromov_witten_degree(Z, [1], 0) == 0
  @test gromov_witten_invariant(X, [pt, pt], [1]; bundle=zero_bundle(X)) == 1
  @test_throws ArgumentError gromov_witten_invariant(X, [pt], [1]; bundle=E)
  @test_throws ArgumentError expected_gromov_witten_degree(
    X, [1], 2; bundle=line_bundle(X, -1)
  )
  @test_throws ArgumentError expected_gromov_witten_degree(
    X, [1], 2; bundle=3 * line_bundle(X, 1)
  )
  @test_throws ArgumentError expected_gromov_witten_degree(X, [1], 2; bundle="O(1)")
  @test_throws ArgumentError expected_gromov_witten_degree(
    X, [1], 2; bundle=structure_sheaf(X)
  )
  @test_throws OverflowError expected_gromov_witten_degree(X, [typemax(Int)], 2)
  @test_throws OverflowError GWF._vertex_bound([typemax(Int)])
  P1 = projective_space(1)
  @test gromov_witten_invariant(P1, [], [1]) == 1
  @test gromov_witten_invariant(P1, [point_class(P1)], [1]) == 1
end
