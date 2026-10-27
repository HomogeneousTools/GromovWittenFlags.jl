"""Recover parent-to-child edges from a preorder depth sequence whose root has depth one."""
function _tree_edges(levels::Vector{Int})
  edges = NTuple{2,Int}[]
  for v in 2:length(levels)
    parent = findlast(i -> i < v && levels[i] == levels[v] - 1, eachindex(levels))
    push!(edges, (parent, v))
  end
  edges
end

"""Build adjacency lists of `(neighbor, edge_index)` pairs, keeping edge labels stable."""
function _tree_adjacency(n::Int, edges::Vector{NTuple{2,Int}})
  adjacency = [Tuple{Int,Int}[] for _ in 1:n]
  for (i, (u, v)) in enumerate(edges)
    push!(adjacency[u], (v, i))
    push!(adjacency[v], (u, i))
  end
  adjacency
end

"""A tree automorphism and the inverse permutation it induces on edge labels."""
struct _TreeAutomorphism
  vertices::Vector{Int}
  edge_inverse::Vector{Int}
end

# Large trees use recursive canonicalization rather than factorial enumeration.
"""Read-only small-tree permutation groups keyed by depth sequence; `nothing` selects the fallback."""
const _TREE_AUTOMORPHISM_CACHE = Dict{
  Tuple{Vararg{Int}},Union{Nothing,Vector{_TreeAutomorphism}}
}()
"""Protect all structural tree caches; reentrant because shape construction requests automorphisms."""
const _TREE_AUTOMORPHISM_LOCK = ReentrantLock()

"""
    _tree_automorphisms(levels, edges)

Cache all uncolored automorphisms for trees with at most eight vertices.
Return `nothing` above this cutoff, where recursive canonicalization avoids
factorial work. Stored permutations are shared read-only by graph callbacks.
"""
function _tree_automorphisms(levels::Vector{Int}, edges::Vector{NTuple{2,Int}})
  key = Tuple(levels)
  lock(_TREE_AUTOMORPHISM_LOCK) do
    get!(_TREE_AUTOMORPHISM_CACHE, key) do
      n = length(levels)
      n > 8 && return nothing
      edge_lookup = Dict(
        minmax(u, v) => edge for (edge, (u, v)) in enumerate(edges)
      )
      result = _TreeAutomorphism[]
      for permutation in permutations(1:n)
        edge_inverse = zeros(Int, length(edges))
        valid = true
        for (edge, (u, v)) in enumerate(edges)
          image = get(edge_lookup, minmax(permutation[u], permutation[v]), 0)
          if iszero(image)
            valid = false
            break
          end
          edge_inverse[image] = edge
        end
        valid && push!(result, _TreeAutomorphism(collect(permutation), edge_inverse))
      end
      result
    end
  end
end

"""Restrict an uncolored automorphism group to those preserving colors and an optional root mark."""
function _color_automorphisms(automorphisms, colors; root::Int=0)
  filter(automorphisms) do automorphism
    (iszero(root) || automorphism.vertices[root] == root) &&
      all(colors[v] == colors[automorphism.vertices[v]] for v in eachindex(colors))
  end
end

"""
    _degree_orbit(degrees, automorphisms)

Return the lexicographically least degree tuple and its stabilizer order.
`automorphisms` is nonempty and must already preserve vertex colors and any
distinguished root. The inverse edge permutation pulls each degree assignment back to the
original edge order; a tuple key owns its values independently of the DFS buffer.
"""
function _degree_orbit(degrees, automorphisms)
  best = first(automorphisms).edge_inverse
  stabilizer = 0
  for automorphism in automorphisms
    inverse = automorphism.edge_inverse
    preserves = true
    better = false
    decided = false
    for edge in eachindex(degrees)
      candidate = degrees[inverse[edge]]
      candidate == degrees[edge] || (preserves = false)
      if !decided && candidate != degrees[best[edge]]
        better = candidate < degrees[best[edge]]
        decided = true
      end
    end
    preserves && (stabilizer += 1)
    better && (best = inverse)
  end
  Tuple(degrees[edge] for edge in best), stabilizer
end

"""Find the one or two centers of a nonempty tree by peeling leaves layer by layer."""
function _tree_centers(adjacency)
  n = length(adjacency)
  n <= 2 && return collect(1:n)
  degrees = length.(adjacency)
  leaves = findall(==(1), degrees)
  remaining = n
  while remaining > 2
    new_leaves = Int[]
    for v in leaves
      degrees[v] = 0
      for (w, _) in adjacency[v]
        if degrees[w] > 1
          degrees[w] -= 1
          degrees[w] == 1 && push!(new_leaves, w)
        end
      end
    end
    remaining -= length(leaves)
    leaves = new_leaves
  end
  leaves
end

"""
    _tree_encoding(adjacency, degrees, colors, marks, vertex, parent)

Return `(encoding, automorphism_order)` for the subtree rooted at `vertex`,
excluding the branch through `parent` (zero at a root). The encoding records
edge degrees, vertex colors, and ordered lists of labeled marks. Equal child
encodings can be permuted freely, contributing their multiplicity factorials.
Each `marks[v]` lists mark indices in increasing order.
The order is a `BigInt`; inputs remain unchanged.
"""
function _tree_encoding(adjacency, degrees, colors, marks, v::Int, parent::Int)
  # Subtree encodings have different tuple types at different depths. This
  # heterogeneous storage is confined to the large-tree canonical fallback.
  children = Tuple{Int,Tuple}[]
  automorphisms = big(1)
  for (w, edge) in adjacency[v]
    w == parent && continue
    encoding, child_automorphisms = _tree_encoding(adjacency, degrees, colors, marks, w, v)
    push!(children, (degrees[edge], encoding))
    automorphisms *= child_automorphisms
  end
  sort!(children)
  i = 1
  while i <= length(children)
    j = i + 1
    while j <= length(children) && children[j] == children[i]
      j += 1
    end
    automorphisms *= factorial(big(j - i))
    i = j
  end
  (colors[v], Tuple(marks[v]), Tuple(children)), automorphisms
end

"""
    _canonical_tree(adjacency, degrees, colors, marks)

Return a relabeling-invariant key and the order of the decoration-preserving
automorphism group. A central edge can be flipped only when its two rooted
halves have identical encodings. BigInt factorials keep large symmetry orders
exact in the fallback used beyond the small-tree enumeration cutoff.
"""
function _canonical_tree(adjacency, degrees, colors, marks)
  centers = _tree_centers(adjacency)
  if length(centers) == 1
    root, automorphisms = _tree_encoding(
      adjacency, degrees, colors, marks, only(centers), 0
    )
    return (0, root), automorphisms
  end
  u, v = centers
  left, left_automorphisms = _tree_encoding(adjacency, degrees, colors, marks, u, v)
  right, right_automorphisms = _tree_encoding(adjacency, degrees, colors, marks, v, u)
  edge = only(i for (w, i) in adjacency[u] if w == v)
  automorphisms = left_automorphisms * right_automorphisms
  left == right && (automorphisms *= 2)
  a, b = minmax(left, right)
  (1, degrees[edge], a, b), automorphisms
end
