"""
    _positive_degree_single_pair_blocks(evaluations, classes, pair, beta, anchor_first)

Two-point localization using a marked root when one insertion has small
fixed-point support. The root carries the anchor insertion; the other labeled
mark is summed over all vertices. Root-preserving automorphisms already
account for the first mark, and orbit–stabilizer accounts for the second.
"""
function _positive_degree_single_pair_blocks(evaluations, classes, pair, beta,
  anchor_first::Bool)
  i, j = pair
  anchor, other = anchor_first ? (i, j) : (j, i)
  data = first(evaluations).data
  T = eltype(first(evaluations).point)
  tables = [
    _restriction_table(evaluation, classes[[anchor, other]]) for evaluation in evaluations
  ]
  anchor_colors = findall(!iszero, @view first(tables)[:, 1])
  totals = [_ExactSum(T) for _ in evaluations]
  max_vertices = _vertex_bound(beta)
  mark_weights = [Vector{T}(undef, max_vertices) for _ in evaluations]
  for n in 2:max_vertices, levels in _FreeTreeIterator(n), color in anchor_colors
    _foreach_rooted_decorated_tree(data, levels, beta, color) do tree
      for k in eachindex(evaluations)
        evaluation = evaluations[k]
        table = tables[k]
        vertex_factor = _vertex_product!(mark_weights[k], evaluation, tree)
        other_sum = sum(
          mark_weights[k][v] * table[tree.colors[v], 2]
          for v in eachindex(tree.colors); init=zero(T))
        iszero(other_sum) && continue
        anchor_value = table[tree.colors[1], 1] * mark_weights[k][1]
        _add_exact!(totals[k],
          vertex_factor * _edge_product(evaluation, tree) * anchor_value * other_sum)
      end
    end
  end
  map(accumulator -> [_sum_value(accumulator)], totals)
end

"""
    _positive_degree_pair_blocks(evaluations, classes, pairs, beta)

Evaluate all requested two-point invariants while sharing unmarked graphs and
restriction tables. On a fixed graph each mark contributes an independent
vertex sum. Their product sums every labeled placement; division by the
unmarked automorphism order gives exactly the marked-orbit multiplicities.
Only vertices and insertion pairs with nonzero restrictions reach edge-factor
arithmetic. A singleton pair with small support uses the rooted algorithm.
`classes` maps the block's local indices to full Schubert indices; each pair
indexes that map. Return one vector per evaluation, in `pairs` order. The
input geometry, class map, pairs, and curve coordinates are read-only.
"""
function _positive_degree_pair_blocks(evaluations, classes, pairs, beta)
  data = first(evaluations).data
  if length(pairs) == 1
    i, j = only(pairs)
    first_evaluation = first(evaluations)
    support_i = count(color -> !iszero(_billey(first_evaluation, classes[i], color)),
      eachindex(data.representatives))
    support_j = count(color -> !iszero(_billey(first_evaluation, classes[j], color)),
      eachindex(data.representatives))
    if min(support_i, support_j) <= max(2, length(data.representatives) ÷ 4)
      return _positive_degree_single_pair_blocks(evaluations, classes, only(pairs),
        beta, support_i <= support_j)
    end
  end
  T = eltype(first(evaluations).point)
  totals = [[_ExactSum(T) for _ in pairs] for _ in evaluations]
  # Both slots use one union of insertions. A symmetric pair needs one table
  # column and one vertex sum, instead of duplicate left/right arithmetic.
  distinct = unique(vcat(first.(pairs), last.(pairs)))
  positions = Dict(class => i for (i, class) in enumerate(distinct))
  pair_positions = [(positions[i], positions[j]) for (i, j) in pairs]
  tables = [_restriction_table(evaluation, classes[distinct]) for evaluation in evaluations]
  max_vertices = _vertex_bound(beta)
  mark_weights = [Vector{T}(undef, max_vertices) for _ in evaluations]
  sums = [zeros(T, length(distinct)) for _ in evaluations]
  for n in 2:max_vertices, levels in _FreeTreeIterator(n)
    _foreach_unmarked_decorated_tree(data, levels, beta) do tree
      for k in eachindex(evaluations)
        evaluation = evaluations[k]
        table = tables[k]
        local_sums = sums[k]
        fill!(local_sums, zero(T))
        vertex_factor = _vertex_product!(mark_weights[k], evaluation, tree)
        for j in eachindex(distinct), vertex in eachindex(tree.colors)
          restriction = table[tree.colors[vertex], j]
          iszero(restriction) && continue
          local_sums[j] += restriction * mark_weights[k][vertex]
        end
        any(pair -> !iszero(local_sums[pair[1]]) && !iszero(local_sums[pair[2]]),
          pair_positions) || continue
        base = vertex_factor * _edge_product(evaluation, tree)
        iszero(base) && continue
        for (index, (i, j)) in enumerate(pair_positions)
          restriction = local_sums[i] * local_sums[j]
          iszero(restriction) || _add_exact!(totals[k][index], base * restriction)
        end
      end
    end
  end
  map(accumulators -> map(_sum_value, accumulators), totals)
end

"""
Preflight poles and compare a complete two-point block at two exact FLINT torus
points. Return a matrix of `Rational{BigInt}` coefficients for quantum assembly.
"""
function _checked_pair_block(data::_GWData, classes, pairs, beta, weights)
  values = _checked_torus_evaluations(data, weights, QQFieldElem) do evaluations
    foreach(evaluation -> _prepare_edges!(evaluation, beta), evaluations)
    _positive_degree_pair_blocks(evaluations, classes, pairs, beta)
  end
  # Keep FLINT rationals inside localization. Matrix assembly, public results,
  # and exact linear algebra retain their Rational{BigInt} coefficient type.
  _BigRational.(values)
end
