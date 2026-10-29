# Retain the original, unpruned coloring iterator as a reference. Graph keys
# and stabilizer orders come from AHU canonicalization, independently of the
# production coloring-prefix comparisons and small-tree permutation orbits.
function decorated_signatures(data, levels, beta, root; pruned)
  found = Any[]
  consume = function (tree)
    marks = [i == 1 && !iszero(root) ? [1] : Int[] for i in eachindex(tree.colors)]
    key, aut = GWF._canonical_tree(tree.adjacency, tree.degrees, tree.colors, marks)
    @test aut == tree.automorphisms
    push!(found, (key, aut))
  end
  if pruned
    if iszero(root)
      GWF._foreach_unmarked_decorated_tree(consume, data, levels, beta)
    else
      GWF._foreach_rooted_decorated_tree(consume, data, levels, beta, root)
    end
  else
    shapes = iszero(root) ? [levels] :
             GWF._rooted_level_sequences(GWF._tree_shape(levels))
    for rooted_levels in shapes
      shape = GWF._tree_shape(rooted_levels)
      tree = GWF._tree_workspace(shape)
      minimum_class = zeros(Int, length(beta))
      iterator = if iszero(root)
        GWF._ColoringIterator(shape.levels, data.neighbors)
      else
        GWF._ColoringIterator(shape.levels, data.neighbors;
          fixed_root=true, root_color=root)
      end
      for colors in iterator
        GWF._color_tree!(minimum_class, tree, data, colors)
        GWF._foreach_colored_degree(consume, tree, shape, data, beta, minimum_class;
          root=iszero(root) ? 0 : 1)
      end
    end
  end
  @test length(found) == length(unique(found))
  Set(found)
end

@testset "Curve-budget pruning preserves decorated orbits" begin
  examples = (projective_space(2), Gr(2, 4), full_flag_variety(TypeA{2}),
    product(projective_space(1), projective_space(2)))
  for X in examples
    data = GWF._gw_data(X)
    for coordinates in Iterators.product(ntuple(_ -> 0:3, picard_rank(X))...)
      1 <= sum(coordinates) <= 3 || continue
      beta = collect(coordinates)
      for n in 2:(1 + sum(beta)), levels in GWF._FreeTreeIterator(n),
        root in [0; collect(1:min(3, length(data.representatives)))]

        @test decorated_signatures(data, levels, beta, root; pruned=true) ==
          decorated_signatures(data, levels, beta, root; pruned=false)
      end
      @test beta == collect(coordinates)
    end
  end
end

@testset "Six-vertex multigraded decorated orbits" begin
  examples = (product(ntuple(_ -> projective_space(1), 4)...),
    product(full_flag_variety(TypeA{2}), projective_space(1), projective_space(1)))
  for X in examples
    data = GWF._gw_data(X)
    for coordinates in Iterators.product(ntuple(_ -> 0:5, 4)...)
      sum(coordinates) == 5 || continue
      beta = collect(coordinates)
      for n in 2:6, levels in GWF._FreeTreeIterator(n), root in (0, 1, 2)
        @test decorated_signatures(data, levels, beta, root; pruned=true) ==
          decorated_signatures(data, levels, beta, root; pruned=false)
      end
    end
  end
end
