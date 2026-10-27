# Canonical graph-homomorphism iterator adapted from GKMtools.jl.  It emits
# one coloring of a free tree by the GKM graph in each isomorphism class,
# instead of generating all root assignments and discarding almost all of
# them after AHU canonicalization.
"""
    _ColoringIterator(levels, neighbors; fixed_root=false, root_color=1)

Enumerate graph homomorphisms from a preorder tree to a GKM graph, modulo
uncolored tree automorphisms. Neighbor lists must be sorted, nonempty for
nontrivial trees, and contain no self loops. Every returned color vector is
owned by its consumer; successor state is kept separate.

Identical sibling subtrees are ordered lexicographically. An isomorphism
exchanging the two central halves imposes a further root-color ordering.
`fixed_root=true` disables this central flip and preserves vertex one's mark.
The reverse traversal and subtree bounds guide which suffix can change next.
"""
struct _ColoringIterator
  levels::Vector{Int}
  neighbors::Vector{Vector{Int}}
  reverse_dfs::Vector{Int}
  reverse_positions::Vector{Int}
  parents::Vector{Int}
  has_central_involution::Bool
  subtree_ends_reverse::Vector{Int}
  subtree_ends::Vector{Int}
  left_siblings::Vector{Int}
  fixed_root::Bool
  root_color::Int
end

"""Last preorder index in the subtree rooted at `vertex`, including the root itself."""
function _subtree_end(levels::Vector{Int}, vertex::Int)
  for i in (vertex + 1):length(levels)
    levels[i] <= levels[vertex] && return i - 1
  end
  length(levels)
end

"""Boundary vertex used to backtrack through a subtree in reverse depth-first order."""
function _subtree_end_reverse(levels, reverse_dfs, vertex)
  position = findfirst(==(vertex), reverse_dfs)
  found = findfirst(
    index -> levels[index] <= levels[vertex],
    @view(reverse_dfs[(position + 1):end]),
  )
  found === nothing ? reverse_dfs[end] : reverse_dfs[position + found]
end

"""Immediate preceding sibling root, or zero when this vertex has no left sibling."""
function _left_sibling_root(parents, vertex)
  parent = parents[vertex]
  iszero(parent) && return 0
  something(findlast(i -> parents[i] == parent, 1:(vertex - 1)), 0)
end

"""Whether the depth sequence admits an involution exchanging identical central halves."""
function _has_central_involution(levels)
  levels == [1, 2] && return true
  isodd(length(levels)) && return false
  second_branch = findfirst(i -> i > 2 && levels[i] == 2, eachindex(levels))
  second_branch === nothing && return false
  @view(levels[3:(second_branch - 1)]) == 1 .+ @view(levels[second_branch:end])
end

function _ColoringIterator(
  levels::Vector{Int}, neighbors::Vector{Vector{Int}};
  fixed_root::Bool=false, root_color::Int=1,
)
  n = length(levels)
  parents = zeros(Int, n)
  for vertex in 2:n
    parents[vertex] = findlast(
      i -> i < vertex && levels[i] == levels[vertex] - 1, eachindex(levels)
    )
  end
  stack = [1]
  reverse_dfs = Int[]
  while !isempty(stack)
    vertex = pop!(stack)
    push!(reverse_dfs, vertex)
    append!(stack, findall(==(vertex), parents))
  end
  subtree_ends = [_subtree_end(levels, vertex) for vertex in 1:n]
  _ColoringIterator(
    levels,
    neighbors,
    reverse_dfs,
    invperm(reverse_dfs),
    parents,
    !fixed_root && _has_central_involution(levels),
    [_subtree_end_reverse(levels, reverse_dfs, vertex) for vertex in 1:n],
    subtree_ends,
    [_left_sibling_root(parents, vertex) for vertex in 1:n],
    fixed_root,
    root_color,
  )
end

"""Lexicographically first valid coloring of a rooted subtree with prescribed root color."""
function _first_coloring(
  levels::Vector{Int}, neighbors::Vector{Vector{Int}}, root_color::Int=1
)
  colors = fill(root_color, length(levels))
  for vertex in 2:length(levels)
    parent = findlast(i -> levels[i] == levels[vertex] - 1, 1:(vertex - 1))
    colors[vertex] = first(neighbors[colors[parent]])
  end
  colors
end

Base.eltype(::Type{_ColoringIterator}) = Vector{Int}
Base.IteratorSize(::Type{_ColoringIterator}) = Base.SizeUnknown()

function Base.iterate(iterator::_ColoringIterator, state=nothing)
  if state === nothing
    colors = _first_coloring(
      iterator.levels, iterator.neighbors, iterator.root_color
    )
    return copy(colors), colors
  end
  colors = state
  look_end = iterator.subtree_ends[1]
  colorable = 0
  position = 2
  while position <= length(iterator.levels)
    vertex = iterator.reverse_dfs[position]
    vertex_end = iterator.subtree_ends[vertex]
    reverse_end = iterator.subtree_ends_reverse[vertex]
    sibling = iterator.left_siblings[vertex]
    sibling_range = iszero(sibling) ? (1:0) : (sibling:(vertex - 1))
    if @view(iterator.levels[vertex:vertex_end]) == @view(iterator.levels[sibling_range]) &&
      @view(colors[vertex:vertex_end]) == @view(colors[sibling_range])
      # Equal sibling blocks must advance together to avoid duplicate orbits.
      position = iterator.reverse_positions[sibling]
      continue
    end
    available = iterator.neighbors[colors[iterator.parents[vertex]]]
    if colors[vertex] < last(available)
      colorable = max(colorable, vertex)
      look_end = reverse_end
    end
    vertex == look_end && !iszero(colorable) && break
    position += 1
  end

  if iszero(colorable)
    iterator.fixed_root && return nothing
    new_root = 0
    if iterator.has_central_involution
      # Choose one orientation of the edge joining the identical central
      # halves: the first endpoint must have the smaller color.
      for candidate in (colors[1] + 1):length(iterator.neighbors)
        candidate < last(iterator.neighbors[candidate]) || continue
        new_root = candidate
        break
      end
    elseif colors[1] < length(iterator.neighbors)
      new_root = colors[1] + 1
    end
    iszero(new_root) && return nothing
    colors[1] = new_root
    colorable = 1
  else
    available = iterator.neighbors[colors[iterator.parents[colorable]]]
    colors[colorable] = available[findfirst(>(colors[colorable]), available)]
  end

  vertex = colorable + 1
  while vertex <= length(iterator.levels)
    root_color = if iterator.has_central_involution && vertex == 2
      available = iterator.neighbors[colors[1]]
      available[findfirst(>(colors[1]), available)]
    else
      first(iterator.neighbors[colors[iterator.parents[vertex]]])
    end
    subtree_end = iterator.subtree_ends[vertex]
    colors[vertex:subtree_end] = _first_coloring(
      iterator.levels[vertex:subtree_end], iterator.neighbors, root_color
    )
    vertex = subtree_end + 1
  end
  copy(colors), colors
end
