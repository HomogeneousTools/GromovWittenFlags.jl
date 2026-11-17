"""
    _BranchGeometry

Degree-truncated GKM geometry for recursive localization. Each edge type is
a root and a covering degree. `valid[a,d]` indexes the remaining-degree
rectangle after attaching edge type `a`; zero means outside that rectangle.
The exponential edge series is shared by the two orientations of an edge.
Cauchy kernels depend only on the distinct incident inverse tangent weights.
"""
struct _BranchGeometry{T}
  degrees::Matrix{Int}
  totals::Vector{Int}
  order::Vector{Int}
  splits::Vector{Vector{NTuple{2,Int}}}
  positive_splits::Vector{Vector{NTuple{2,Int}}}
  roots::Vector{Int}
  covers::Vector{Int}
  costs::Matrix{Int}
  sizes::Vector{Int}
  valid::Matrix{Int}
  remainder::Matrix{Int}
  endpoints::Matrix{Int}
  reverse_types::Matrix{Int}
  groups::Matrix{Int}
  sources::Vector{Int}
  targets::Vector{Int}
  group_types::Vector{Int}
  inverse_weights::Matrix{T}
  edge_weights::Matrix{T}
  weight_ids::Matrix{Int}
  distinct_weights::Vector{T}
  kernels::Matrix{T}
  kernel_known::BitMatrix
end

"""Add a product to an owned accumulator; cached factors are never mutated."""
_branch_add_product!(total, a, b) = total + a * b
_branch_add_product!(total::QQFieldElem, a::QQFieldElem, b::QQFieldElem) =
  addmul!(total, a, b)

"""Return mixed-radix strides, checking dimension products before allocating arrays."""
function _branch_strides(bounds)
  strides = Int[]
  size = 1
  for b in bounds
    push!(strides, size)
    size = Base.Checked.checked_mul(size, Base.Checked.checked_add(b, 1))
  end
  strides, size
end

"""Return a degree's one-based mixed-radix index."""
_branch_degree_index(degree, strides) =
  1 + sum(degree[i] * strides[i] for i in eachindex(strides))

"""Obtain an incident-weight Cauchy kernel, checking all poles before recursion."""
function _branch_kernel!(geometry::_BranchGeometry{T}, i::Int, j::Int) where {T}
  if !geometry.kernel_known[i, j]
    u, v = geometry.distinct_weights[i], geometry.distinct_weights[j]
    reciprocal = _divide_or_pole(one(T), u + v)
    geometry.kernels[i, j] = v^2 * reciprocal
    geometry.kernels[j, i] = u^2 * reciprocal
    geometry.kernel_known[i, j] = geometry.kernel_known[j, i] = true
  end
  geometry.kernels[i, j]
end

"""Build edge, degree-convolution, and pole-preflight data without listing any trees."""
function _branch_geometry(
  evaluation::_GWEvaluation{D,T}, beta; progress=nothing
) where {D,T}
  _vertex_bound(beta)
  data = evaluation.data
  strides, n_degrees = _branch_strides(beta)
  p = length(beta)
  degrees = Matrix{Int}(undef, p, n_degrees)
  for (i, degree) in enumerate(CartesianIndices(Tuple(beta .+ 1)))
    degrees[:, i] .= Tuple(degree) .- 1
  end
  totals = vec(sum(degrees; dims=1))
  order = sortperm(1:n_degrees; by=i -> (totals[i], i))
  splits = [NTuple{2,Int}[] for _ in 1:n_degrees]
  positive = [NTuple{2,Int}[] for _ in 1:n_degrees]
  for d in 1:n_degrees
    degree = @view degrees[:, d]
    for gamma in CartesianIndices(Tuple(degree .+ 1))
      left = collect(Tuple(gamma) .- 1)
      pair = (_branch_degree_index(left, strides),
        _branch_degree_index(degree - left, strides))
      push!(splits[d], pair)
      pair[1] == 1 || push!(positive[d], pair)
    end
  end
  roots, covers = Int[], Int[]
  costs_list = Vector{Int}[]
  lookup = Dict{Tuple{Int,Int},Int}()
  for (root, curve) in enumerate(data.curve_classes)
    maximum = minimum(beta[i] ÷ curve[i] for i in eachindex(beta) if curve[i] > 0)
    for cover in 1:maximum
      push!(roots, root)
      push!(covers, cover)
      push!(costs_list, cover .* curve)
      lookup[(root, cover)] = length(roots)
    end
  end
  n_types = length(roots)
  costs = isempty(costs_list) ? zeros(Int, p, 0) : hcat(costs_list...)
  sizes = zeros(Int, n_types)
  valid, remainder = zeros(Int, n_types, n_degrees), zeros(Int, n_types, n_degrees)
  for a in 1:n_types
    local_strides, sizes[a] = _branch_strides(beta - costs[:, a])
    for d in 1:n_degrees
      degree = @view degrees[:, d]
      if all(degree[i] + costs[i, a] <= beta[i] for i in 1:p)
        valid[a, d] = _branch_degree_index(degree, local_strides)
      end
      if all(degree[i] >= costs[i, a] for i in 1:p)
        remainder[a, d] = _branch_degree_index(degree - costs[:, a], strides)
      end
    end
  end
  n_colors = length(data.representatives)
  endpoints, reverse_types = zeros(Int, n_colors, n_types), zeros(Int, n_colors, n_types)
  inverse_weights = Matrix{T}(undef, n_colors, n_types)
  weight_ids = zeros(Int, n_colors, n_types)
  distinct_weights = T[]
  distinct_ids = Dict{_BigRational,Int}()
  for color in 1:n_colors, a in 1:n_types
    root, cover = roots[a], covers[a]
    endpoint = data.steps[color, root]
    endpoints[color, a] = endpoint
    reverse_types[color, a] = lookup[(data.edge_roots[endpoint][color], cover)]
    alpha = evaluation.root_values[color, data.nonparabolic_indices[root]]
    u = _divide_or_pole(-T(cover), alpha)
    inverse_weights[color, a] = u
    uid = get!(distinct_ids, _BigRational(u)) do
      push!(distinct_weights, u)
      length(distinct_weights)
    end
    weight_ids[color, a] = uid
  end
  groups = zeros(Int, n_colors, n_types)
  sources, targets, group_types = Int[], Int[], Int[]
  edge_weights = Matrix{T}(undef, n_colors, n_types)
  for color in 1:n_colors, a in 1:n_types
    groups[color, a] == 0 || continue
    target, reverse = endpoints[color, a], reverse_types[color, a]
    @assert costs[:, a] == costs[:, reverse]
    @assert inverse_weights[color, a] == -inverse_weights[target, reverse]
    push!(sources, color)
    push!(targets, target)
    push!(group_types, a)
    groups[color, a] = groups[target, reverse] = length(sources)
    value = _edge_factor(evaluation, color, roots[a], covers[a])
    edge_weights[color, a] = edge_weights[target, reverse] = value
  end
  n_weights = length(distinct_weights)
  geometry = _BranchGeometry(degrees, totals, order, splits, positive, roots, covers,
    costs, sizes, valid, remainder, endpoints, reverse_types, groups, sources, targets,
    group_types, inverse_weights, edge_weights, weight_ids, distinct_weights,
    Matrix{T}(undef, n_weights, n_weights), falses(n_weights, n_weights))
  # A compatible pair can occur at a bivalent vertex. Detect numerical poles
  # now, rather than losing a long recursion to a late exceptional torus point.
  for a in 1:n_types, b in a:n_types
    all(costs[i, a] + costs[i, b] <= beta[i] for i in 1:p) || continue
    for color in 1:n_colors
      _branch_kernel!(geometry, weight_ids[color, a], weight_ids[color, b])
    end
  end
  progress === nothing || progress((phase=:geometry, degree_states=n_degrees,
    fixed_points=n_colors, oriented_edge_types=n_types,
    distinct_inverse_weights=n_weights))
  geometry
end

"""Disjoint label-subset convolutions for nilpotent marking variables."""
function _branch_mask_splits(n_marks::Int)
  0 <= n_marks <= 7 ||
    throw(ArgumentError("recursive localization supports at most eight labeled marks"))
  result = [NTuple{2,Int}[] for _ in 1:(1 << n_marks)]
  for mask in 0:(length(result) - 1)
    part = mask
    while true
      push!(result[mask + 1], (part + 1, xor(mask, part) + 1))
      part == 0 && break
      part = (part - 1) & mask
    end
  end
  result
end

"""
    _branch_root_values(evaluation, mark_values, beta; geometry, progress)

Sum all positive-degree trees rooted at a distinguished marked vertex.
`mark_values` gives the restrictions of the other labeled insertions. The
returned vector is paired with the distinguished insertion's restrictions.

For incident inverse weights u, normalized vertex Euler values e, and edge
weights K, the complete branch sum obeys the formal series identities

  tau_i = t_i + e_i * sum_a u_a*C_ia,
  D_i(u) = sum_a u_a^2*C_ia/(u+u_a),
  Z_a = exp(u_a*(tau_i-tau_j)),
  C_ia = K_a*Q^degree(a)*Z_a*(-1/u_a-e_j*D_j(-u_a)).

The root contribution is sum_a C_ia + (t_i*tau_i-tau_i^2/2)/e_i.
Lagrange inversion of the vertex exponential generating function proves
these identities, including leaf and unstable bivalent contributions.
Every C coefficient depends only on strictly smaller degrees; Z is computed
after tau. Equal edge orientations share Z. No free trees are enumerated.

All cached coefficients are read-only. Mutating FLINT accumulation is used
only on fresh local temporaries, never on shared zeros or stored factors.
"""
function _branch_root_values(evaluation::_GWEvaluation{D,T}, mark_values,
  beta; geometry=nothing, progress=nothing) where {D,T}
  geometry === nothing && (geometry = _branch_geometry(evaluation, beta; progress))
  g = geometry
  n_colors, n_marks = size(mark_values)
  n_colors == length(evaluation.data.representatives) ||
    throw(DimensionMismatch("mark restrictions"))
  mask_splits = _branch_mask_splits(n_marks)
  n_masks = length(mask_splits)
  n_degrees = length(g.totals)
  n_types = length(g.roots)
  # Filling with a shared zero is safe: stored entries are replaced, never mutated.
  branch = [fill(zero(T), n_masks, g.sizes[a], n_colors) for a in 1:n_types]
  exponential = [fill(zero(T), n_masks, g.sizes[a]) for a in g.group_types]
  tau = fill(zero(T), n_masks, n_degrees, n_colors)
  for color in 1:n_colors, mark in 1:n_marks
    tau[(1 << (mark - 1)) + 1, 1, color] = mark_values[color, mark]
  end
  for edge in eachindex(g.sources)
    source, target, a = g.sources[edge], g.targets[edge], g.group_types[edge]
    u = g.inverse_weights[source, a]
    exponential[edge][1, 1] = one(T)
    for mask in 1:(n_masks - 1)
      bit = trailing_zeros(mask)
      previous = xor(mask, 1 << bit)
      exponential[edge][mask + 1, 1] =
        exponential[edge][previous + 1, 1] *
        u * (mark_values[source, bit + 1] - mark_values[target, bit + 1])
    end
  end
  current = fill(zero(T), n_masks, n_types)
  root_values = zeros(T, n_colors)
  final_degree = n_degrees # The upper corner is last in mixed-radix order.
  started = time_ns()
  for d in g.order
    d == 1 && continue
    active = findall(a -> g.remainder[a, d] != 0, 1:n_types)
    incoming = findall(a -> g.valid[a, d] != 0, 1:n_types)
    for color in 1:n_colors
      fill!(current, zero(T))
      euler = evaluation.vertex_euler[color]
      for a in active
        target, reverse = g.endpoints[color, a], g.reverse_types[color, a]
        edge = g.groups[color, a]
        eta = g.remainder[a, d]
        u = g.inverse_weights[color, a]
        factor = g.edge_weights[color, a]
        for mask in 1:n_masks
          total = zero(T)
          for (left, right) in g.positive_splits[eta]
            li, ri = g.valid[reverse, left], g.valid[a, right]
            for (lm, rm) in mask_splits[mask]
              x, y = branch[reverse][lm, li, target], exponential[edge][rm, ri]
              iszero(x) || iszero(y) || (total = _branch_add_product!(total, x, y))
            end
          end
          current[mask, a] =
            factor *
            (
              -exponential[edge][mask, g.valid[a, eta]] / u -
              evaluation.vertex_euler[target] * total
            )
        end
      end
      for mask in 1:n_masks
        total = zero(T)
        for a in active
          value = current[mask, a]
          iszero(value) ||
            (total = _branch_add_product!(total, g.inverse_weights[color, a], value))
        end
        tau[mask, d, color] = euler * total
      end
      # Project the current outgoing messages to all needed incoming states.
      for a in incoming
        li = g.valid[a, d]
        for mask in 1:n_masks
          total = zero(T)
          for b in active
            value = current[mask, b]
            iszero(value) && continue
            kernel = g.kernels[g.weight_ids[color, a], g.weight_ids[color, b]]
            total = _branch_add_product!(total, kernel, value)
          end
          branch[a][mask, li, color] = total
        end
      end
      if d == final_degree
        mask = n_masks
        root = sum((current[mask, a] for a in active); init=zero(T))
        local_marks = zero(T)
        full_mask = mask - 1
        for mark in 1:n_marks
          reduced = xor(full_mask, 1 << (mark - 1)) + 1
          local_marks = _branch_add_product!(local_marks,
            mark_values[color, mark], tau[reduced, d, color])
        end
        square = zero(T)
        for (left, right) in g.splits[d], (lm, rm) in mask_splits[mask]
          x, y = tau[lm, left, color], tau[rm, right, color]
          iszero(x) || iszero(y) || (square = _branch_add_product!(square, x, y))
        end
        root_values[color] = root + (local_marks - square / 2) / euler
      end
    end
    # Euler differentiation computes exp(u*(tau_i-tau_j)) by convolution.
    for edge in eachindex(g.sources)
      source, target, a = g.sources[edge], g.targets[edge], g.group_types[edge]
      li = g.valid[a, d]
      li == 0 && continue
      u = g.inverse_weights[source, a]
      for mask in 1:n_masks
        total = zero(T)
        for (left, right) in g.positive_splits[d]
          ri = g.valid[a, right]
          for (lm, rm) in mask_splits[mask]
            difference = tau[lm, left, source] - tau[lm, left, target]
            value = exponential[edge][rm, ri]
            iszero(difference) || iszero(value) ||
              (total = _branch_add_product!(total, g.totals[left] * difference, value))
          end
        end
        exponential[edge][mask, li] = u * total / g.totals[d]
      end
    end
    progress === nothing || progress((phase=:degree, degree=collect(g.degrees[:, d]),
      degree_index=d, degree_states=n_degrees, mark_channels=n_masks,
      seconds=(time_ns() - started) / 1e9))
  end
  progress === nothing || progress((phase=:complete, degree_states=n_degrees,
    mark_channels=n_masks,
    stored_coefficients=n_masks *
                        (
      n_colors * sum(g.sizes) + sum(g.sizes[a] for a in g.group_types) +
      n_colors * n_degrees
    ), seconds=(time_ns() - started) / 1e9))
  root_values
end

"""Compute one positive-degree primary invariant by rooting its first insertion."""
function _branch_invariant(
  evaluation::_GWEvaluation{D,T}, classes, beta; progress=nothing
) where {D,T}
  if isempty(classes)
    # The divisor equation recovers an unmarked invariant using an ample divisor.
    data = evaluation.data
    divisor_indices = [
      _class_index(data, data.generators[node]) for node in marked_nodes(data.variety)
    ]
    root = [
      sum((_billey(evaluation, index, color) for index in divisor_indices); init=zero(T))
      for color in eachindex(data.representatives)
    ]
    values = _branch_root_values(
      evaluation, Matrix{T}(undef, length(root), 0), beta; progress
    )
    return sum(root .* values) / T(sum(beta))
  end
  table = _restriction_table(evaluation, classes)
  values = _branch_root_values(evaluation, table[:, 2:end], beta; progress)
  sum(table[:, 1] .* values)
end

"""Choose a small set of second-insertion sources covering a two-point block."""
function _branch_pair_sources(pairs)
  sources = Set(i for (i, j) in pairs if i == j)
  remaining = [(i, j) for (i, j) in pairs if !(i in sources || j in sources)]
  while !isempty(remaining)
    choices = unique(vcat(first.(remaining), last.(remaining)))
    source = first(sort!(choices; by=i -> (-count(pair -> i in pair, remaining), i)))
    push!(sources, source)
    filter!(pair -> !(pair[1] in sources || pair[2] in sources), remaining)
  end
  sort!(collect(sources))
end

"""Sum a two-point block by branch recursion, sharing geometry between sources."""
function _branch_pair_blocks(evaluation::_GWEvaluation{D,T}, classes, pairs, beta;
  progress=nothing) where {D,T}
  geometry = _branch_geometry(evaluation, beta; progress)
  table = _restriction_table(evaluation, classes)
  sources = _branch_pair_sources(pairs)
  assignments = [j in sources ? j : i for (i, j) in pairs]
  values = zeros(T, length(pairs))
  for source in sources
    callback =
      progress === nothing ? nothing : event -> progress(merge(event, (source=source,)))
    roots = _branch_root_values(evaluation, table[:, source:source], beta;
      geometry, progress=callback)
    for (index, (i, j)) in enumerate(pairs)
      assignments[index] == source || continue
      other = i == source ? j : i
      values[index] = sum(table[:, other] .* roots)
    end
  end
  values
end
