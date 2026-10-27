"""
    _TreeShape

Read-only combinatorics of a preorder tree: edges, incidence, a parent-before-
child traversal, and its small-tree automorphism group. Shapes are independent
of the flag variety, curve class, torus point, and defining bundle.
"""
struct _TreeShape
  levels::Vector{Int}
  edges::Vector{NTuple{2,Int}}
  adjacency::Vector{Vector{Tuple{Int,Int}}}
  parent_vertex::Vector{Int}
  parent_edge::Vector{Int}
  order::Vector{Int}
  automorphisms::Union{Nothing,Vector{_TreeAutomorphism}}
end

"""Shared read-only tree incidence and traversal data, keyed by owned depth-sequence tuples."""
const _TREE_SHAPE_CACHE = Dict{Tuple{Vararg{Int}},_TreeShape}()
"""Rooted depth sequences, one per vertex orbit; protected by the structural tree lock."""
const _ROOTED_LEVEL_CACHE = Dict{Tuple{Vararg{Int}},Vector{Vector{Int}}}()

"""
    _DecoratedTree

Borrowed callback workspace for an invariant-map fixed graph. Edges are
oriented from their parent in `order`; each `roots[edge]` is expressed in the
parent color's root frame, and `degrees[edge]` is its positive cover degree.
The automorphism order preserves colors, degrees, and a root mark when used.

Enumeration reuses the mutable buffers. A callback must copy any data it
retains beyond its call. There is no producer task or Channel to leak when
localization throws at an exceptional torus point.
"""
mutable struct _DecoratedTree
  edges::Vector{NTuple{2,Int}}
  adjacency::Vector{Vector{Tuple{Int,Int}}}
  roots::Vector{Int}
  degrees::Vector{Int}
  colors::Vector{Int}
  parent_vertex::Vector{Int}
  parent_edge::Vector{Int}
  order::Vector{Int}
  automorphisms::Union{Int,BigInt}
end

"""Orient a tree from `root` and return parent vertices, parent edges, and traversal order."""
function _rooted_tree(adjacency, root::Int=1)
  n = length(adjacency)
  parent_vertex = zeros(Int, n)
  parent_edge = zeros(Int, n)
  order = [root]
  for v in order, (w, edge) in adjacency[v]
    w == parent_vertex[v] && continue
    parent_vertex[w] = v
    parent_edge[w] = edge
    push!(order, w)
  end
  parent_vertex, parent_edge, order
end

"""Return a shared tree shape, copying the iterator's depth sequence into the cache."""
function _tree_shape(levels::Vector{Int})
  lock(_TREE_AUTOMORPHISM_LOCK) do
    get!(_TREE_SHAPE_CACHE, Tuple(levels)) do
      owned = copy(levels)
      edges = _tree_edges(owned)
      adjacency = _tree_adjacency(length(owned), edges)
      parents, parent_edges, order = _rooted_tree(adjacency)
      _TreeShape(owned, edges, adjacency, parents, parent_edges, order,
        _tree_automorphisms(owned, edges))
    end
  end
end

"""Allocate the mutable color/root/degree buffers borrowed by one shape's graph callbacks."""
function _tree_workspace(shape::_TreeShape)
  _DecoratedTree(shape.edges, shape.adjacency, zeros(Int, length(shape.edges)),
    ones(Int, length(shape.edges)), zeros(Int, length(shape.levels)),
    shape.parent_vertex, shape.parent_edge, shape.order, 1)
end

"""
    _foreach_degree_assignment(f, data, roots, beta)

Stream positive cover degrees satisfying `sum(d[e] * curve[roots[e]]) == beta`.
The callback borrows one degree buffer. Suffix sums reserve at least one cover
for each remaining edge, pruning impossible branches before recursion. The
last edge is solved directly rather than trying every possible degree.
"""
function _foreach_degree_assignment(f::F, data::_GWData, roots, beta) where {F}
  n = length(roots)
  current = zeros(Int, n)
  suffix = zeros(Int, length(beta), n + 1)
  for edge in n:-1:1, i in eachindex(beta)
    suffix[i, edge] = suffix[i, edge + 1] + data.curve_classes[roots[edge]][i]
  end
  remaining = copy(beta)
  function visit(edge::Int)
    if edge > n
      all(iszero, remaining) && f(current)
      return nothing
    end
    curve = data.curve_classes[roots[edge]]
    maximum_degree = typemax(Int)
    for i in eachindex(beta)
      budget = remaining[i] - suffix[i, edge + 1]
      budget < 0 && return nothing
      curve[i] > 0 && (maximum_degree = min(maximum_degree, budget ÷ curve[i]))
    end
    maximum_degree == typemax(Int) && return nothing
    if edge == n
      maximum_degree > 0 || return nothing
      all(remaining[i] == maximum_degree * curve[i] for i in eachindex(beta)) ||
        return nothing
      current[edge] = maximum_degree
      f(current)
      return nothing
    end
    for degree in 1:maximum_degree
      current[edge] = degree
      for i in eachindex(remaining)
        remaining[i] -= degree * curve[i]
      end
      visit(edge + 1)
      for i in eachindex(remaining)
        remaining[i] += degree * curve[i]
      end
    end
    nothing
  end
  visit(1)
  nothing
end

"""Fill edge roots from a graph coloring and return its minimal (all covers degree one) class."""
function _color_tree!(minimum_class, tree::_DecoratedTree, data::_GWData, colors)
  tree.colors = colors
  fill!(minimum_class, 0)
  for vertex in @view tree.order[2:end]
    edge = tree.parent_edge[vertex]
    parent_color = colors[tree.parent_vertex[vertex]]
    root = data.edge_roots[parent_color][colors[vertex]]
    tree.roots[edge] = root
    curve = data.curve_classes[root]
    for i in eachindex(minimum_class)
      minimum_class[i] += curve[i]
    end
  end
  nothing
end

"""
    _foreach_colored_degree(f, tree, shape, data, beta, minimum_class; root=0)

Decorate a fixed coloring by degree orbits and compute each stabilizer order.
Small shapes use their cached permutation group; larger shapes use recursive
canonical encodings. A distinguished root restricts automorphisms to those
fixing that vertex. The all-unit case avoids allocating any orbit keys.
"""
function _foreach_colored_degree(f::F, tree::_DecoratedTree, shape::_TreeShape,
  data::_GWData, beta, minimum_class; root::Int=0) where {F}
  any(i -> minimum_class[i] > beta[i], eachindex(beta)) && return nothing
  preserving = if shape.automorphisms === nothing
    nothing
  else
    _color_automorphisms(shape.automorphisms, tree.colors; root)
  end
  marks = preserving === nothing ? [Int[] for _ in tree.colors] : nothing
  if !iszero(root) && marks !== nothing
    push!(marks[root], 1)
  end
  if minimum_class == beta
    fill!(tree.degrees, 1)
    tree.automorphisms = if preserving === nothing
      last(_canonical_tree(tree.adjacency, tree.degrees, tree.colors, marks))
    else
      length(preserving)
    end
    f(tree)
    return nothing
  end
  seen = Set{Tuple}()
  _foreach_degree_assignment(data, tree.roots, beta) do degrees
    key, aut = if preserving === nothing
      _canonical_tree(tree.adjacency, degrees, tree.colors, marks)
    else
      _degree_orbit(degrees, preserving)
    end
    key in seen && return nothing
    push!(seen, key)
    tree.degrees = degrees
    tree.automorphisms = aut
    f(tree)
    nothing
  end
  nothing
end

"""
    _foreach_unmarked_decorated_tree(f, data, levels, beta)

Stream every colored, degree-decorated free-tree orbit in curve class `beta`,
without marks. Canonical colorings remove uncolored symmetries first; degree
orbits then quotient only by the color stabilizer. The callback borrows all
workspace arrays, and exceptions propagate synchronously to the caller.
"""
function _foreach_unmarked_decorated_tree(f::F, data::_GWData, levels,
  beta::Vector{Int}) where {F}
  shape = _tree_shape(levels)
  tree = _tree_workspace(shape)
  minimum_class = zeros(Int, length(beta))
  for colors in _ColoringIterator(shape.levels, data.neighbors)
    _color_tree!(minimum_class, tree, data, colors)
    _foreach_colored_degree(f, tree, shape, data, beta, minimum_class)
  end
  nothing
end

"""Encode an uncolored tree rooted at a vertex, sorting child shapes into a canonical preorder."""
function _rooted_encoding(adjacency, vertex::Int, parent::Int)
  children = [
    _rooted_encoding(adjacency, child, vertex)
    for (child, _) in adjacency[vertex] if child != parent
  ]
  sort!(children; by=first)
  key = Tuple(first(child) for child in children)
  levels = Int[1]
  for (_, child_levels) in children
    append!(levels, child_levels .+ 1)
  end
  key, levels
end

"""Cache one rooted depth sequence for each vertex orbit of an uncolored free tree."""
function _rooted_level_sequences(shape::_TreeShape)
  lock(_TREE_AUTOMORPHISM_LOCK) do
    get!(_ROOTED_LEVEL_CACHE, Tuple(shape.levels)) do
      sequences = Dict{Tuple,Vector{Int}}()
      for root in eachindex(shape.levels)
        key, levels = _rooted_encoding(shape.adjacency, root, 0)
        get!(sequences, key, levels)
      end
      sort!(collect(values(sequences)))
    end
  end
end

"""
    _foreach_rooted_decorated_tree(f, data, levels, beta, root_color)

Stream rooted decorated-tree orbits whose distinguished vertex has the given
fixed-point color. Each uncolored vertex orbit supplies one rooted shape;
color and degree stabilizers must preserve the root mark. The distinguished
vertex is always vertex one in the normalized preorder.
"""
function _foreach_rooted_decorated_tree(f::F, data::_GWData, levels,
  beta::Vector{Int}, root_color::Int) where {F}
  for rooted_levels in _rooted_level_sequences(_tree_shape(levels))
    shape = _tree_shape(rooted_levels)
    tree = _tree_workspace(shape)
    minimum_class = zeros(Int, length(beta))
    for colors in
        _ColoringIterator(shape.levels, data.neighbors; fixed_root=true, root_color)
      _color_tree!(minimum_class, tree, data, colors)
      _foreach_colored_degree(f, tree, shape, data, beta, minimum_class; root=1)
    end
  end
  nothing
end
