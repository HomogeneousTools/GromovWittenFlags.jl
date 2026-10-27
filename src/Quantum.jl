"""
    _classical_matrices(data, dim, chern_coordinates, weights)

Compute the full ambient metric and first-Chern three-point pairing together.
At each of two torus points, Schubert restrictions and the Euler density are
computed once and reused for every entry. Codimension gates exclude classes
that cannot contribute on the target. The returned matrices are symmetric.
"""
function _classical_matrices(data::_GWData, dim::Int, c, weights)
  active = findall(degree -> degree <= dim, data.codimensions)
  metric_pairs = [
    (i, j) for i in eachindex(active) for j in i:length(active)
    if data.codimensions[active[i]] + data.codimensions[active[j]] == dim
  ]
  chern_pairs = [
    (i, j) for i in eachindex(active) for j in i:length(active)
    if data.codimensions[active[i]] + data.codimensions[active[j]] == dim - 1
  ]
  divisors = [
    (c[node], _class_index(data, data.generators[node]))
    for node in marked_nodes(data.variety) if !iszero(c[node])
  ]
  _checked_torus_evaluations(data, weights) do evaluations
    map(evaluations) do evaluation
      T = eltype(evaluation.point)
      n = length(data.representatives)
      metric = zeros(T, n, n)
      chern = zeros(T, n, n)
      table = _restriction_table(evaluation, active)
      for color in 1:n
        density = inv(evaluation.vertex_euler[color])
        chern_value = sum(
          coefficient * _billey(evaluation, index, color)
          for (coefficient, index) in divisors; init=zero(T))
        for (pairs, matrix, multiplier) in ((metric_pairs, metric, density),
          (chern_pairs, chern, density * chern_value))
          iszero(multiplier) && continue
          for (i, j) in pairs
            value = table[color, i] * table[color, j]
            iszero(value) && continue
            a, b = active[i], active[j]
            matrix[a, b] += multiplier * value
            a == b || (matrix[b, a] = matrix[a, b])
          end
        end
      end
      metric, chern
    end
  end
end

"""
    _compact_beta_vectors(coefficients, bound)

Enumerate nonnegative integer vectors of weighted degree at most `bound`.
All coefficients must be strictly positive, ensuring a finite list. Results
are ordered by total coordinate sum, then lexicographically, with zero first.
"""
function _compact_beta_vectors(coefficients::Vector{Int}, bound::Int)
  result = Vector{Vector{Int}}()
  current = zeros(Int, length(coefficients))
  function visit(i::Int, remaining::Int)
    if i > length(coefficients)
      push!(result, copy(current))
      return nothing
    end
    for value in 0:(remaining ÷ coefficients[i])
      current[i] = value
      visit(i + 1, remaining - value * coefficients[i])
    end
    nothing
  end
  visit(1, bound)
  sort!(result; by=beta -> (sum(beta), beta))
end

"""Enumerate ambient Fano degrees from validated target data without expanding bundle weights again."""
function _fano_curve_classes(data::_GWData, dim::Int, c)
  coefficients = Int[c[node] for node in marked_nodes(data.variety)]
  any(<(0), coefficients) && throw(
    ArgumentError(
      "automatic curve enumeration requires nef first Chern class; supply explicit betas"
    ),
  )
  index = foldl(gcd, coefficients; init=0)
  all(iszero, coefficients) && return index, [zeros(Int, length(coefficients))]
  any(iszero, coefficients) && throw(
    ArgumentError(
      "automatic curve enumeration is unbounded in a zero Chern direction; supply explicit betas"
    ),
  )
  index, _compact_beta_vectors(coefficients, dim + 1)
end

"""
    fano_curve_classes(X; bundle=nothing)
    fano_curve_classes(Z::ZeroLocus)

Return `(index, betas)` for the flag variety or homogeneous zero locus.
For strictly positive first Chern coordinates, enumerate every nonnegative
ambient curve class with first Chern degree at most `dimension + 1`, including
zero. The index is the gcd of the ambient first Chern coordinates.
If all coordinates vanish, return only zero. Mixed zero and positive, or
negative, coordinates require explicit `betas` in
[`small_quantum_multiplication`](@ref).
"""
function fano_curve_classes(X::PartialFlagVariety; bundle=nothing)
  data = _gw_data(X)
  dim, c, _ = _complete_intersection_data(data, bundle)
  _fano_curve_classes(data, dim, c)
end

"""
    small_quantum_multiplication(X; bundle=nothing, betas=nothing)
    small_quantum_multiplication(Z::ZeroLocus; betas=nothing)

Compute the flag-ambient projection of small quantum multiplication by the
target's first Chern class. Return `(matrix, grading, basis, metric)`.
Matrix entries are exact [`NovikovPolynomial`](@ref)s; [`evaluate`](@ref)
specializes their variables. Columns are input classes and rows are output
coefficients. All returned basis words and polynomials are owned by the caller.

For a zero locus, ambient Schubert restrictions are reduced modulo the radical
of the intersection pairing. The ambient subspace need not be preserved by the
full quantum product; nonambient cohomology is excluded. Grading is complex
codimension minus half the target dimension.

Explicit `betas` define a finite truncation: include zero for the classical
term; duplicates are ignored. Otherwise [`fano_curve_classes`](@ref) supplies
the finite list. The defining bundle must satisfy the supported Euler-twist
weight condition described in the manual.
"""
function small_quantum_multiplication(X::PartialFlagVariety; bundle=nothing, betas=nothing)
  _small_quantum_multiplication(_gw_data(X), bundle, betas)
end

"""Validate target and curve inputs before constructing intersection matrices."""
function _small_quantum_multiplication(data::_GWData, bundle, betas)
  dim, c, weights = _complete_intersection_data(data, bundle)
  beta_list = if betas === nothing
    last(_fano_curve_classes(data, dim, c))
  else
    unique([_normalize_beta(data, beta) for beta in betas])
  end
  _quantum_multiplication(Val(picard_rank(data.variety)), data, dim, c, weights, beta_list)
end

"""
    _quantum_multiplication(Val(P), data, dim, c, weights, beta_list)

Assemble multiplication degree by degree. One exact LU factorization raises
an index of each rational pairing matrix; sparse polynomial terms are added
only after this solve. This avoids a matrix inverse and repeated polynomial
copies in a cubic matrix-multiplication loop. `P` specializes the Novikov rank.
"""
function _quantum_multiplication(::Val{P}, data::_GWData, dim::Int, c, weights,
  beta_list) where {P}
  metric_full, classical_full = _classical_matrices(data, dim, c, weights)
  indices = if weights === nothing
    collect(eachindex(data.representatives))
  else
    _independent_indices(metric_full)
  end
  metric = metric_full[indices, indices]
  classical = classical_full[indices, indices]
  basis = [_owned_schubert(data, data.representatives[i]) for i in indices]
  factorization = lu(metric)
  n = length(indices)
  matrix = [zero(NovikovPolynomial{P}) for _ in 1:n, _ in 1:n]
  for beta in beta_list
    coefficients = if all(iszero, beta)
      classical
    else
      c_beta = _curve_degree(data, c, beta)
      iszero(c_beta) && continue
      required_degree = Base.Checked.checked_add(dim - 1, c_beta)
      pairs = [
        (i, j) for i in 1:n for j in i:n
        if
        data.codimensions[indices[i]] + data.codimensions[indices[j]] == required_degree
      ]
      isempty(pairs) && continue
      invariants = _checked_pair_block(data, indices, pairs, beta, weights)
      block = zeros(_BigRational, n, n)
      for ((i, j), invariant) in zip(pairs, invariants)
        block[i, j] = c_beta * invariant
        block[j, i] = block[i, j]
      end
      block
    end
    values = factorization \ coefficients
    for column in 1:n, row in 1:n
      iszero(values[row, column]) && continue
      _add_novikov_term!(matrix[row, column], beta, values[row, column])
    end
  end
  grading = Diagonal([(2 * data.codimensions[i] - dim) // 2 for i in indices])
  (; matrix, grading, basis, metric)
end

fano_curve_classes(Z::ZeroLocus) =
  fano_curve_classes(ambient_variety(Z); bundle=defining_bundle(Z))

function small_quantum_multiplication(Z::ZeroLocus; betas=nothing)
  small_quantum_multiplication(ambient_variety(Z); bundle=defining_bundle(Z), betas)
end
