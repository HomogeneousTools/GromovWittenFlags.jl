@testset "Tree enumeration and symmetry factors" begin
  @test_throws ArgumentError GWF._FreeTreeIterator(0)
  @test [length(collect(GWF._FreeTreeIterator(n))) for n in 1:8] ==
    [1, 1, 1, 2, 3, 6, 11, 23]
  for n in 1:8
    trees = collect(GWF._FreeTreeIterator(n))
    @test length(unique(Tuple.(trees))) == length(trees)
    @test all(levels -> length(levels) == n && first(levels) == 1, trees)
  end
  # Nine vertices cross the explicit-permutation cutoff. A star's symmetry
  # order is independently known: every leaf permutation is an automorphism.
  levels = [1; fill(2, 8)]
  shape = GWF._tree_shape(levels)
  @test shape.automorphisms === nothing
  marks = [Int[] for _ in levels]
  colors = [1; fill(2, 8)]
  @test last(GWF._canonical_tree(shape.adjacency, ones(Int, 8), colors, marks)) ==
    factorial(big(8))
  data = GWF._gw_data(projective_space(1))
  orders = Int[]
  GWF._foreach_unmarked_decorated_tree(data, levels, [8]) do tree
    push!(orders, tree.automorphisms)
  end
  @test orders == fill(factorial(8), 2)

  # The decorated fallback must distinguish a double cover, reducing the leaf
  # stabilizer from 8! to 7!, without counting all eight placements separately.
  orders = Int[]
  GWF._foreach_unmarked_decorated_tree(data, levels, [9]) do tree
    push!(orders, tree.automorphisms)
  end
  @test orders == fill(factorial(7), 2)
  for n in 2:6, levels in GWF._FreeTreeIterator(n)
    edges = GWF._tree_edges(levels)
    adjacency = GWF._tree_adjacency(n, edges)
    colors = [mod(i, 2) + 1 for i in 1:n]
    degrees = [mod(i, 2) + 1 for i in eachindex(edges)]
    marks = [Int[] for _ in 1:n]
    key, aut = GWF._canonical_tree(adjacency, degrees, colors, marks)
    automorphisms = GWF._tree_automorphisms(levels, edges)
    brute = count(automorphisms) do a
      all(colors[i] == colors[a.vertices[i]] for i in 1:n) &&
        all(degrees[i] == degrees[a.edge_inverse[i]] for i in eachindex(edges))
    end
    @test aut == brute
    preserving = GWF._color_automorphisms(automorphisms, colors)
    degree_key, stabilizer = GWF._degree_orbit(degrees, preserving)
    @test degree_key == minimum(Tuple(degrees[a.edge_inverse]) for a in preserving)
    @test stabilizer == brute
    # Relabel vertices in reverse order; canonicalization must be invariant.
    relabeled = [(n + 1 - u, n + 1 - v) for (u, v) in edges]
    other = GWF._tree_adjacency(n, relabeled)
    @test GWF._canonical_tree(other, degrees, reverse(colors), reverse(marks)) == (key, aut)
  end

  graphs = ([[2, 3], [1, 3], [1, 2]], [[2], [1, 3], [2, 4], [3]])
  for neighbors in graphs, n in 2:6, levels in GWF._FreeTreeIterator(n)
    edges = GWF._tree_edges(levels)
    adjacency = GWF._tree_adjacency(n, edges)
    degrees = ones(Int, n - 1)
    marks = [Int[] for _ in 1:n]
    brute = Set{Any}()
    for assignment in Iterators.product(ntuple(_ -> eachindex(neighbors), n)...)
      all(assignment[v] in neighbors[assignment[u]] for (u, v) in edges) || continue
      push!(
        brute, first(GWF._canonical_tree(adjacency, degrees, collect(assignment), marks))
      )
    end
    colorings = collect(GWF._ColoringIterator(levels, neighbors))
    canonical = [
      first(GWF._canonical_tree(adjacency, degrees, c, marks)) for c in colorings
    ]
    @test Set(canonical) == brute
    @test length(unique(canonical)) == length(canonical)
    @test length(unique(objectid.(colorings))) == length(colorings)
  end
end

@testset "Rooted coloring orbits and owned outputs" begin
  graphs = ([[2, 3], [1, 3], [1, 2]], [[2], [1, 3], [2, 4], [3]])
  for neighbors in graphs, n in 1:5, unrooted in GWF._FreeTreeIterator(n),
    levels in GWF._rooted_level_sequences(GWF._tree_shape(unrooted))

    edges = GWF._tree_edges(levels)
    adjacency = GWF._tree_adjacency(n, edges)
    degrees = ones(Int, n - 1)
    # Distinguish vertex one so canonicalization allows only root-preserving
    # isomorphisms, independently of the iterator's sibling ordering.
    marks = [i == 1 ? [1] : Int[] for i in 1:n]
    for root_color in eachindex(neighbors)
      brute = Set{Any}()
      for assignment in Iterators.product(ntuple(_ -> eachindex(neighbors), n)...)
        first(assignment) == root_color || continue
        all(assignment[v] in neighbors[assignment[u]] for (u, v) in edges) || continue
        push!(brute,
          first(GWF._canonical_tree(adjacency, degrees, collect(assignment), marks)))
      end
      iterator = GWF._ColoringIterator(levels, neighbors; fixed_root=true, root_color)
      colorings = collect(iterator)
      canonical = [
        first(GWF._canonical_tree(adjacency, degrees, colors, marks))
        for colors in colorings
      ]
      @test Set(canonical) == brute
      @test length(unique(canonical)) == length(canonical)
    end
  end

  iterator = GWF._ColoringIterator([1, 2, 2], first(graphs))
  expected = collect(iterator)
  colors, state = iterate(iterator)
  @test colors !== state
  fill!(colors, 0)
  successor, _ = iterate(iterator, state)
  @test successor == expected[2]
end
