# Literal reduced-subword enumeration is deliberately separate from the
# production automaton: it makes no use of weak-order states or transitions.
function subword_restriction(evaluation, class_index, color)
  data = evaluation.data
  target = data.representatives[class_index]
  letters = GWF.word(data.representatives[color])
  total = zero(QQ)
  for positions in combinations(eachindex(letters), length(target))
    element = one(data.weyl_group)
    value = one(QQ)
    for position in positions
      element *= data.generators[Int(letters[position])]
      value *= evaluation.word_values[color][position]
    end
    element == target && (total += value)
  end
  total
end

# Brute-force labeled mark placements, followed by canonicalization, provide
# an independent check on the orbit–stabilizer summation used in production.
function marked_orbit_invariant(evaluation, classes, beta)
  total = zero(QQ)
  for n in 2:(sum(beta) + 1), levels in GWF._FreeTreeIterator(n)
    GWF._foreach_unmarked_decorated_tree(evaluation.data, levels, beta) do tree
      mark_weights = zeros(QQ, n)
      base = GWF._flag_tree_unmarked_factor!(mark_weights, evaluation, tree)
      seen = Set{Tuple}()
      for positions in Iterators.product(ntuple(_ -> 1:n, length(classes))...)
        marks = [Int[] for _ in 1:n]
        for (mark, vertex) in enumerate(positions)
          push!(marks[vertex], mark)
        end
        key, aut = GWF._canonical_tree(tree.adjacency, tree.degrees, tree.colors, marks)
        key in seen && continue
        push!(seen, key)
        contribution = base * tree.automorphisms / aut
        for (mark, vertex) in enumerate(positions)
          contribution *=
            mark_weights[vertex] *
            GWF._billey(evaluation, classes[mark], tree.colors[vertex])
        end
        total += contribution
      end
    end
  end
  total
end

@testset "Schubert recurrence versus literal subwords" begin
  for X in (full_flag_variety(TypeA{3}), full_flag_variety(TypeB{2}),
    full_flag_variety(TypeG2), Gr(2, 4), product(projective_space(1), projective_space(1)))
    data = GWF._gw_data(X)
    evaluation = GWF._gw_evaluation(data, 0)
    for color in eachindex(data.representatives), class in eachindex(data.representatives)
      @test GWF._billey(evaluation, class, color) ==
        subword_restriction(evaluation, class, color)
    end
    for color in eachindex(data.representatives), root in eachindex(data.curve_classes)
      endpoint = data.steps[color, root]
      reverse_root = data.edge_roots[endpoint][color]
      @test evaluation.root_values[color, data.nonparabolic_indices[root]] ==
        -evaluation.root_values[endpoint, data.nonparabolic_indices[reverse_root]]
    end
  end
end

@testset "Labeled marks and vertex cache" begin
  P1 = projective_space(1)
  P2 = projective_space(2)
  F = full_flag_variety(TypeA{2})
  P3 = projective_space(3)
  cases = (
    (P1, [], [1], nothing),
    (P1, [point_class(P1)], [1], nothing),
    (P1, fill(point_class(P1), 2), [1], nothing),
    (P1, fill(point_class(P1), 2), [2], nothing),
    (P2, [schubert_class(P2, [1]), point_class(P2), point_class(P2)], [1], nothing),
    (P2, fill(point_class(P2), 5), [2], nothing),
    (F, fill(point_class(F), 2), [1, 1], nothing),
    (P3, [], [1], line_bundle(P3, 3)),
  )
  for (X, insertions, beta, bundle) in cases
    data = GWF._gw_data(X)
    _, _, weights = GWF._complete_intersection_data(data, bundle)
    evaluation = GWF._gw_evaluation(data, 0; weights)
    GWF._prepare_edges!(evaluation, beta)
    classes = [GWF._class_index(data, class) for class in insertions]
    @test only(GWF._positive_degree_invariants((evaluation,), classes, beta)) ==
      marked_orbit_invariant(evaluation, classes, beta)
    for n in 2:(sum(beta) + 1), levels in GWF._FreeTreeIterator(n)
      GWF._foreach_unmarked_decorated_tree(data, levels, beta) do tree
        for vertex in eachindex(tree.colors)
          @test GWF._vertex_factors(evaluation, tree, vertex) ==
            GWF._vertex_factors_uncached(evaluation, tree, vertex)
        end
      end
    end
  end

  data = GWF._gw_data(P1)
  evaluation = GWF._gw_evaluation(data, 0)
  tree = GWF._tree_workspace(GWF._tree_shape([1, 2]))
  tree.colors .= [1, 2]
  tree.roots .= [1]
  @test GWF._flag_tree_unmarked_factor!(zeros(QQ, 2), evaluation, tree) == 1
  tree.degrees .= [2]
  alpha = evaluation.root_values[1, data.nonparabolic_indices[1]]
  @test GWF._flag_tree_unmarked_factor!(zeros(QQ, 2), evaluation, tree) ==
    -1 / (2 * alpha^2)
end

@testset "Exact poles and retry policy" begin
  data = GWF._gw_data(projective_space(1))
  @test_throws GWF._LocalizationPole GWF._divide_or_pole(one(QQ), zero(QQ))
  @test GWF._divide_or_pole(big(3)//big(2), big(5)//big(7)) == 21//10
  calls = Ref(0)
  result = GWF._checked_torus_evaluations(data, nothing) do _
    calls[] += 1
    calls[] == 1 && GWF._divide_or_pole(one(QQ), zero(QQ))
    (big(7)//big(1), big(7)//big(1))
  end
  @test result == 7
  @test calls[] == 2
  calls[] = 0
  @test_throws DivideError GWF._checked_torus_evaluations(data, nothing) do _
    calls[] += 1
    div(1, 0)
  end
  @test calls[] == 1
  calls[] = 0
  @test_throws DomainError GWF._checked_torus_evaluations(data, nothing) do _
    calls[] += 1
    throw(DomainError("programming errors must propagate"))
  end
  @test calls[] == 1
  calls[] = 0
  @test_throws ErrorException GWF._checked_torus_evaluations(data, nothing) do _
    calls[] += 1
    (one(QQ), big(2)//big(1))
  end
  @test calls[] == 1
  calls[] = 0
  result = GWF._checked_torus_evaluations(data, nothing) do _
    calls[] += 1
    calls[] == 1 ? (big(1)//big(0), big(1)//big(0)) : (one(QQ), one(QQ))
  end
  @test result == 1
  @test calls[] == 2
  calls[] = 0
  @test_throws ErrorException GWF._checked_torus_evaluations(data, nothing) do _
    calls[] += 1
    GWF._divide_or_pole(one(QQ), zero(QQ))
  end
  @test calls[] == 8
  @test_throws ErrorException GWF._foreach_unmarked_decorated_tree(data, [1, 2], [1]) do _
    error("callback exceptions propagate without leaving a producer task")
  end
end

@testset "Degree DFS versus Cartesian enumeration" begin
  data = GWF._gw_data(full_flag_variety(TypeA{2}))
  for roots in ([1, 1], [1, 2], [1, 3], [3, 3]),
    beta in ([0, 0], [1, 1], [2, 1], [2, 2], [3, 2])

    actual = Set{Tuple}()
    GWF._foreach_degree_assignment(data, roots, beta) do degrees
      push!(actual, Tuple(degrees))
    end
    expected = Set{Tuple}()
    for degrees in Iterators.product(ntuple(_ -> 1:sum(beta), length(roots))...)
      all(
        sum(degrees[e] * data.curve_classes[roots[e]][i] for e in eachindex(roots)) ==
        beta[i]
        for i in eachindex(beta)
      ) && push!(expected, degrees)
    end
    @test actual == expected
  end
  # The callback sees a reusable buffer. Retaining copies must preserve all
  # assignments, and enumeration must leave caller-supplied inputs alone.
  root = only(findall(==([1, 0]), data.curve_classes))
  roots, beta = [root, root], [4, 0]
  retained = Vector{Int}[]
  GWF._foreach_degree_assignment(data, roots, beta) do degrees
    push!(retained, copy(degrees))
  end
  @test Set(Tuple.(retained)) == Set(((1, 3), (2, 2), (3, 1)))
  @test roots == [root, root] && beta == [4, 0]
  @test length(unique(objectid.(retained))) == 3
end
