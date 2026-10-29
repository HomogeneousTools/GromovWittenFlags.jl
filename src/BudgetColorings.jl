"""
    _budget_color_comparisons(levels, parents)

Precompute identical-sibling prefix comparisons for a normalized preorder
tree. Entry `v` contains `(left, right, length)` comparisons that become
decidable when vertex `v` receives its color. The input arrays are read-only.
"""
function _budget_color_comparisons(levels, parents)
  comparisons = [NTuple{3,Int}[] for _ in levels]
  for right in 2:length(levels)
    left = _left_sibling_root(parents, right)
    iszero(left) && continue
    last = _subtree_end(levels, right)
    length_left = right - left
    last - right + 1 == length_left || continue
    @view(levels[left:(right - 1)]) == @view(levels[right:last]) || continue
    for vertex in right:last
      push!(comparisons[vertex], (left, right, vertex - right + 1))
    end
  end
  comparisons
end

"""Whether an assigned color prefix respects nonincreasing identical-sibling colors."""
function _budget_color_canonical(colors, comparisons, vertex)
  for (left, right, n) in comparisons[vertex]
    for k in 0:(n - 1)
      colors[left + k] > colors[right + k] && break
      colors[left + k] < colors[right + k] && return false
    end
  end
  true
end

"""
    _visit_budget_colors!(f, tree, data, beta, remaining, minimum_class,
                         neighbors, comparisons, central_flip, vertex)

Assign preorder colors starting at `vertex`, subtracting each unit edge class
from `remaining` before descending. All cover degrees are positive, so a
negative coordinate rules out every completion of that prefix. Restore the
budget on normal return; colors, roots, and `minimum_class` are scratch storage.
Call `f(tree, minimum_class)` for each admissible canonical coloring. Both
arguments are borrowed until the callback returns; exceptions propagate.
"""
function _visit_budget_colors!(f::F, tree, data, beta, remaining, minimum_class,
  neighbors, comparisons, central_flip, vertex::Int) where {F}
  if vertex > length(tree.colors)
    for i in eachindex(beta)
      minimum_class[i] = beta[i] - remaining[i]
    end
    f(tree, minimum_class)
    return nothing
  end
  parent = tree.parent_vertex[vertex]
  parent_color = tree.colors[parent]
  edge = tree.parent_edge[vertex]
  for color in neighbors[parent_color]
    central_flip && vertex == 2 && color < parent_color && continue
    root = data.edge_roots[parent_color][color]
    curve = data.curve_classes[root]
    all(curve[i] <= remaining[i] for i in eachindex(beta)) || continue
    tree.colors[vertex] = color
    _budget_color_canonical(tree.colors, comparisons, vertex) || continue
    tree.roots[edge] = root
    for i in eachindex(beta)
      remaining[i] -= curve[i]
    end
    _visit_budget_colors!(f, tree, data, beta, remaining, minimum_class,
      neighbors, comparisons, central_flip, vertex + 1)
    for i in eachindex(beta)
      remaining[i] += curve[i]
    end
  end
  nothing
end

"""
    _foreach_budget_colors(f, tree, shape, data, beta; root_color=0)

Visit canonical colorings whose all-unit edge class fits `beta`. Identical
sibling prefixes and the unmarked central flip use the same conventions as
`_ColoringIterator`. A nonzero `root_color` fixes vertex one and disables the
central flip. The callback borrows `tree` and the minimum class buffer; degree
assignment and stabilizer computation remain the caller's responsibility.
The supplied workspace is overwritten; `shape`, `data`, and `beta` are read-only.
"""
function _foreach_budget_colors(f::F, tree, shape, data, beta;
  root_color::Int=0) where {F}
  comparisons = _budget_color_comparisons(shape.levels, tree.parent_vertex)
  neighbors = [
    [
      v for v in data.neighbors[color]
      if all(
        data.curve_classes[data.edge_roots[color][v]][i] <= beta[i]
        for i in eachindex(beta)
      )
    ]
    for color in eachindex(data.neighbors)
  ]
  remaining = copy(beta)
  minimum_class = zeros(Int, length(beta))
  central_flip = iszero(root_color) && _has_central_involution(shape.levels)
  roots = iszero(root_color) ? eachindex(data.neighbors) : (root_color:root_color)
  for color in roots
    tree.colors[1] = color
    _visit_budget_colors!(f, tree, data, beta, remaining, minimum_class,
      neighbors, comparisons, central_flip, 2)
  end
  nothing
end
