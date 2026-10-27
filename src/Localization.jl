"""
    _GWEvaluation{D,T,W}

One exact torus specialization for one target. Root values and Schubert
restrictions use scalar type `T`; `weights` binds this evaluation to a single
defining bundle, so the Euler-factor caches cannot be reused for another twist.
All caches are local to the computation and require no global arithmetic state.

`vertex_euler[color]` is the tangent Euler class divided by the bundle Euler
class. Edge factors exclude endpoint Euler classes, which have been collected
at vertices according to their valence. This reduces large intermediate
numerators and denominators before summing localization contributions.
"""
struct _GWEvaluation{D,T,W}
  data::D
  point::Vector{T}
  root_values::Matrix{T}
  word_values::Vector{Vector{T}}
  rfactor::Vector{T}
  weights::W
  fiber_values::Matrix{T}
  splittings::Matrix{Int}
  vertex_euler::Vector{T}
  billey::Dict{Tuple{Int,Int},T}
  edge_factor::Dict{NTuple{3,Int},T}
  inverse_flags::Dict{NTuple{3,Int},NTuple{2,T}}
  vertex_factor::Dict{NTuple{5,Int},NTuple{2,T}}
end

"""
    _divide_or_pole(numerator, denominator)

Divide exact scalars, treating a zero denominator as a localization pole.
Julia rational arithmetic permits `1//0`; an explicit check is therefore
necessary before division, rather than relying on a thrown arithmetic error.
"""
function _divide_or_pole(numerator, denominator)
  iszero(denominator) && throw(DivideError())
  numerator / denominator
end

"""Evaluate an integral root at small integral torus coordinates, checking overflow."""
function _root_value(coefficients, point)
  value = 0
  for i in eachindex(point)
    value = Base.Checked.checked_add(value,
      Base.Checked.checked_mul(coefficients[i], point[i]))
  end
  value
end

"""
    _gw_evaluation(data, seed, T=Rational{BigInt}; weights=nothing)

Construct a deterministic torus specialization and its target-specific caches.
A xorshift sequence avoids the weight hyperplanes met by arithmetic
progressions. Weyl actions and reduced-word roots are already cached in `data`;
their integral linear forms are evaluated before conversion to exact rationals.
"""
function _gw_evaluation(data::_GWData, seed::Int, ::Type{T}=_BigRational;
  weights=nothing) where {T}
  R = length(data.generators)
  state = UInt64(0x9e3779b97f4a7c15) * UInt64(seed + 1)
  coordinates = Vector{Int}(undef, R)
  for i in eachindex(coordinates)
    state = xor(state, state << 13)
    state = xor(state, state >> 7)
    state = xor(state, state << 17)
    coordinates[i] = 1009 + Int(state % UInt64(1_000_003))
  end
  point = T.(coordinates)
  values = map(c -> T(_root_value(c, coordinates)), data.root_coordinates)
  word_values = [T[_root_value(c, coordinates) for c in roots] for roots in data.word_roots]
  rfactor = T[
    prod((-values[color, root] for root in data.nonparabolic_indices); init=one(T))
    for color in eachindex(data.representatives)
  ]
  any(iszero, rfactor) && throw(DivideError())

  n_weights = weights === nothing ? 0 : length(weights)
  fiber_values = Matrix{T}(undef, length(rfactor), n_weights)
  splittings = Matrix{Int}(undef, n_weights, length(data.nonparabolic_indices))
  vertex_euler = copy(rfactor)
  if n_weights > 0
    # Fundamental-weight coordinates must be converted to simple-root
    # coordinates before evaluation. The minus sign converts the package's
    # Borel–Weil–Bott weights into actual bundle-fiber characters.
    Cinv = _BigRational.(cartan_matrix_inverse(dynkin_type(data.variety)))
    for (j, weight) in enumerate(weights)
      for (color, w_inverse) in enumerate(data.inverse_representatives)
        root_coordinates = Cinv * coefficients(weight * w_inverse)
        fiber_values[color, j] = -sum(T(root_coordinates[i]) * point[i] for i in 1:R)
      end
      for (q, root) in enumerate(data.nonparabolic_indices)
        splittings[j, q] = _root_coroot_pairing(weight, data.positive_coroots[root])
      end
    end
    for color in eachindex(rfactor)
      bundle_euler = prod(@view fiber_values[color, :])
      vertex_euler[color] = _divide_or_pole(rfactor[color], bundle_euler)
    end
  end
  _GWEvaluation(data, point, values, word_values, rfactor, weights, fiber_values,
    splittings, vertex_euler, Dict{Tuple{Int,Int},T}(), Dict{NTuple{3,Int},T}(),
    Dict{NTuple{3,Int},NTuple{2,T}}(), Dict{NTuple{5,Int},NTuple{2,T}}())
end

"""Resolve an insertion to a cached coset index, rejecting incompatible Dynkin types."""
function _class_index(data::_GWData{DT,R}, class::WeylGroupElem{DT,R}) where {DT,R}
  data.representative_index[coefficients(data.stabilizer_weight * inv(class))]
end

_class_index(::_GWData, class) =
  throw(ArgumentError("insertions must be Weyl group elements of the ambient Dynkin type"))

"""
    _billey(evaluation, class_index, vertex_index)

Evaluate Billey's reduced-subword formula for a Schubert class at a fixed
point. The cached automaton retains only prefixes that can reach the target.
States are updated in descending length, so a reduced-word letter is used
at most once. The diagonal restriction is simply the product of word roots.
"""
function _billey(evaluation::_GWEvaluation{D,T}, class_index::Int,
  vertex_index::Int) where {D,T}
  get!(evaluation.billey, (class_index, vertex_index)) do
    data = evaluation.data
    class_index == 1 && return one(T)
    data.codimensions[vertex_index] < data.codimensions[class_index] && return zero(T)
    betas = evaluation.word_values[vertex_index]
    class_index == vertex_index && return prod(betas; init=one(T))
    data.codimensions[vertex_index] == data.codimensions[class_index] && return zero(T)
    plan = _billey_plan(data, class_index)
    values = zeros(T, size(plan.steps, 1))
    values[1] = one(T)
    for (s, beta) in zip(word(data.representatives[vertex_index]), betas)
      for i in (length(values) - 1):-1:1
        successor = plan.steps[i, Int(s)]
        iszero(successor) && continue
        iszero(values[i]) && continue
        values[successor] += beta * values[i]
      end
    end
    last(values)
  end
end

"""
    _edge_factor(evaluation, color, root, degree)

Moving-edge localization contribution, including the cover factor `1/degree`
and the bundle Euler twist, but excluding endpoint tangent Euler classes.
The latter combine at each vertex as `vertex_euler^(valence - 1)`.
The root–coroot pairing in the edge formula is cached before specialization.
"""
function _edge_factor(evaluation::_GWEvaluation{D,T}, color::Int,
  root::Int, degree::Int) where {D,T}
  get!(evaluation.edge_factor, (color, root, degree)) do
    data = evaluation.data
    alpha = evaluation.root_values[color, data.nonparabolic_indices[root]]
    # The tangent direction along the covered curve contributes the factorial
    # quotient. Its cyclic degree-fold deck group supplies the extra 1/degree.
    factor =
      T((-1)^degree) * T(degree)^(2 * degree - 1) /
      (T(factorial(big(degree)))^2 * alpha^(2 * degree))
    for (q, gamma_index) in enumerate(data.nonparabolic_indices)
      q == root && continue
      m = degree * data.root_pairings[q, root]
      gamma = evaluation.root_values[color, gamma_index]
      # Negative splittings contribute moving H¹ weights in the numerator;
      # nonnegative splittings contribute moving H⁰ weights in the denominator.
      if m <= -2
        for j in 1:(-m - 1)
          factor *= -(gamma + (T(j) / T(degree)) * alpha)
        end
      elseif m >= 0
        for k in 0:m
          factor = _divide_or_pole(factor, -(gamma - (T(k) / T(degree)) * alpha))
        end
      end
    end
    # The Euler twist multiplies the section weights, including both endpoints
    # of each nonnegative splitting. Endpoint normalization is done at vertices.
    for j in axes(evaluation.fiber_values, 2)
      value = evaluation.fiber_values[color, j]
      for k in 0:(evaluation.splittings[j, root] * degree)
        factor *= value + (T(k) / T(degree)) * alpha
      end
    end
    factor
  end
end

"""
    _prepare_edges!(evaluation, beta)

Precompute every covered-edge factor that can occur in `beta`. Besides warming
the local cache, this detects edge poles before support pruning can skip a
zero-restriction graph. Every actual edge degree is bounded by its contribution
to each nonzero curve coordinate, so the preflight covers the entire tree sum.
"""
function _prepare_edges!(evaluation::_GWEvaluation, beta)
  _vertex_bound(beta)
  data = evaluation.data
  for (root, curve) in enumerate(data.curve_classes)
    maximum_degree = minimum(beta[i] ÷ curve[i] for i in eachindex(beta) if curve[i] > 0)
    for degree in 1:maximum_degree, color in eachindex(data.representatives)
      _edge_factor(evaluation, color, root, degree)
    end
  end
  nothing
end

"""Cache inverse tangent weights at both ends of a covered invariant curve."""
function _inverse_tangent_weights(evaluation::_GWEvaluation{D,T}, color::Int,
  root::Int, degree::Int) where {D,T}
  get!(evaluation.inverse_flags, (color, root, degree)) do
    alpha = evaluation.root_values[color, evaluation.data.nonparabolic_indices[root]]
    value = _divide_or_pole(T(degree), alpha)
    (-value, value)
  end
end

"""
    _vertex_factors_uncached(evaluation, tree, vertex)

Return `(zero_mark_factor, mark_weight)` at a vertex. If incident tangent
weights are `w_i`, the descendant integral contributes
`prod(1/w_i) * sum(1/w_i)^(valence - 3 + n_marks)`. Each labeled mark therefore
contributes one additional factor `sum(1/w_i)`. The zero-mark factor includes
the normalized tangent/bundle Euler contribution at this vertex.
"""
function _vertex_factors_uncached(evaluation::_GWEvaluation{D,T},
  tree::_DecoratedTree, vertex::Int) where {D,T}
  inverse_product = one(T)
  inverse_sum = zero(T)
  for (_, edge) in tree.adjacency[vertex]
    incoming = tree.parent_edge[vertex] == edge
    parent = incoming ? tree.parent_vertex[vertex] : vertex
    flags = _inverse_tangent_weights(evaluation, tree.colors[parent],
      tree.roots[edge], tree.degrees[edge])
    inverse_weight = flags[incoming ? 2 : 1]
    inverse_product *= inverse_weight
    inverse_sum += inverse_weight
  end
  valence = length(tree.adjacency[vertex])
  zero_marks = if valence < 3
    _divide_or_pole(inverse_product, inverse_sum^(3 - valence))
  else
    inverse_product * inverse_sum^(valence - 3)
  end
  euler = evaluation.vertex_euler[tree.colors[vertex]]
  zero_marks * euler^(valence - 1), inverse_sum
end

"""
    _vertex_factors(evaluation, tree, vertex)

Cache leaf and bivalent-vertex factors by color and the unordered incident
`(local root, cover degree)` pairs. These vertices dominate tree sums and
otherwise repeat the same rational divisions millions of times. Higher
valences are evaluated directly, avoiding a large combinatorial cache.
"""
function _vertex_factors(evaluation::_GWEvaluation, tree::_DecoratedTree, vertex::Int)
  adjacency = tree.adjacency[vertex]
  length(adjacency) > 2 && return _vertex_factors_uncached(evaluation, tree, vertex)
  color = tree.colors[vertex]
  neighbor, edge = adjacency[1]
  first_flag = (
    evaluation.data.edge_roots[color][tree.colors[neighbor]], tree.degrees[edge]
  )
  second_flag = (0, 0)
  if length(adjacency) == 2
    neighbor, edge = adjacency[2]
    second_flag = (
      evaluation.data.edge_roots[color][tree.colors[neighbor]], tree.degrees[edge]
    )
    first_flag, second_flag = minmax(first_flag, second_flag)
  end
  key = (color, first_flag..., second_flag...)
  get!(evaluation.vertex_factor, key) do
    _vertex_factors_uncached(evaluation, tree, vertex)
  end
end

"""Fill borrowed mark weights and multiply all zero-mark vertex contributions."""
function _vertex_product!(mark_weights, evaluation::_GWEvaluation{D,T},
  tree::_DecoratedTree) where {D,T}
  factor = one(T)
  for vertex in eachindex(tree.colors)
    zero_marks, mark_weight = _vertex_factors(evaluation, tree, vertex)
    factor *= zero_marks
    mark_weights[vertex] = mark_weight
  end
  factor
end

"""Multiply cached edge contributions and divide by the decorated-tree automorphism order."""
function _edge_product(evaluation::_GWEvaluation{D,T}, tree::_DecoratedTree) where {D,T}
  factor = one(T) / T(tree.automorphisms)
  for (edge, (parent, _)) in enumerate(tree.edges)
    factor *= _edge_factor(
      evaluation, tree.colors[parent], tree.roots[edge], tree.degrees[edge]
    )
    iszero(factor) && break
  end
  factor
end

"""Complete unmarked graph contribution; `mark_weights` is overwritten up to the vertex count."""
function _flag_tree_unmarked_factor!(
  mark_weights, evaluation::_GWEvaluation, tree::_DecoratedTree
)
  _vertex_product!(mark_weights, evaluation, tree) * _edge_product(evaluation, tree)
end

"""Classical intersection integral at one torus point, including the optional Euler twist."""
function _degree_zero_invariant(evaluation::_GWEvaluation{D,T}, classes) where {D,T}
  total = zero(T)
  for color in eachindex(evaluation.data.representatives)
    term = inv(evaluation.vertex_euler[color])
    for class in classes
      term *= _billey(evaluation, class, color)
      iszero(term) && break
    end
    total += term
  end
  total
end

"""
    _positive_degree_invariants(evaluations, classes, beta)

Compute a primary invariant at several torus points while enumerating each
unmarked graph once. Labeled mark placements factor into independent sums
over vertices. Orbit–stabilizer turns division by marked automorphisms into
division by unmarked automorphisms, eliminating Cartesian mark enumeration.
Repeated insertions share one restriction table and one vertex sum.
"""
function _positive_degree_invariants(evaluations, classes, beta)
  data = first(evaluations).data
  T = eltype(first(evaluations).point)
  distinct = unique(classes)
  multiplicities = [count(==(class), classes) for class in distinct]
  tables = [_restriction_table(evaluation, distinct) for evaluation in evaluations]
  totals = [_ExactSum(T) for _ in evaluations]
  max_vertices = _vertex_bound(beta)
  mark_weights = [Vector{T}(undef, max_vertices) for _ in evaluations]
  for n in 2:max_vertices, levels in _FreeTreeIterator(n)
    _foreach_unmarked_decorated_tree(data, levels, beta) do tree
      for k in eachindex(evaluations)
        evaluation = evaluations[k]
        vertex_factor = _vertex_product!(mark_weights[k], evaluation, tree)
        restriction = one(T)
        for j in eachindex(distinct)
          marked_sum = sum(
            mark_weights[k][v] * tables[k][tree.colors[v], j]
            for v in eachindex(tree.colors); init=zero(T))
          restriction *= marked_sum^multiplicities[j]
          iszero(restriction) && break
        end
        iszero(restriction) && continue
        _add_exact!(
          totals[k], restriction * vertex_factor * _edge_product(evaluation, tree)
        )
      end
    end
  end
  map(_sum_value, totals)
end

"""Test finiteness recursively for scalar, tuple, or matrix-valued localization results."""
_all_finite(value::Number) = isfinite(value)
_all_finite(values) = all(_all_finite, values)

"""
    _checked_torus_evaluations(f, data, weights)

Run `f` on two target-bound exact evaluations and require equal finite results.
Only explicit localization poles (`DivideError`) trigger a retry; unrelated
domain and programming errors propagate. At most eight deterministic pairs
are attempted. A disagreement is a correctness failure, never a retry.
"""
function _checked_torus_evaluations(f::F, data::_GWData, weights) where {F}
  for seed in 0:2:14
    try
      evaluations = (
        _gw_evaluation(data, seed; weights), _gw_evaluation(data, seed + 1; weights)
      )
      values = f(evaluations)
      _all_finite(values) || throw(DivideError())
      values[1] == values[2] ||
        error("localization evaluations disagree at two exact torus points")
      return values[1]
    catch exception
      exception isa DivideError || rethrow()
    end
  end
  error("localization encountered poles at every exact torus-point pair")
end

"""Dispatch degree-zero integration or shared positive-degree graph evaluation with exact checks."""
function _checked_localization(data::_GWData, classes, beta, weights)
  _checked_torus_evaluations(data, weights) do evaluations
    if all(iszero, beta)
      return map(evaluation -> _degree_zero_invariant(evaluation, classes), evaluations)
    end
    foreach(evaluation -> _prepare_edges!(evaluation, beta), evaluations)
    if length(classes) == 2
      blocks = _positive_degree_pair_blocks(evaluations, classes, [(1, 2)], beta)
      map(only, blocks)
    else
      _positive_degree_invariants(evaluations, classes, beta)
    end
  end
end
