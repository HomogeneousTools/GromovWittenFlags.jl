"""
    _classical_matrices(data, dim, chern_coordinates, weights)

Return `(metric, chern_pairing)`, both full ambient symmetric matrices with
Schubert-class indices on both axes. `dim`, `chern_coordinates`, and `weights`
describe the validated target, rather than necessarily its ambient variety.
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
  _append_curve_vectors!(result, current, coefficients, 1, bound)
  sort!(result; by=beta -> (sum(beta), beta))
end

"""
    _append_curve_vectors!(result, current, coefficients, position, remaining)

Append owned curve vectors with weighted degree at most the remaining budget.
The prefix before `position` is fixed; the worker overwrites the remaining
coordinates of `current`. Coefficients must be positive. Copying only complete
vectors preserves result ownership without allocating at each search branch.
"""
function _append_curve_vectors!(
  result, current, coefficients, position::Int, remaining::Int
)
  if position > length(coefficients)
    push!(result, copy(current))
    return nothing
  end
  for value in 0:(remaining ÷ coefficients[position])
    current[position] = value
    _append_curve_vectors!(result, current, coefficients, position + 1,
      remaining - value * coefficients[position])
  end
  nothing
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

Return the tuple `(index::Int, betas::Vector{Vector{Int}})` for the flag variety
or homogeneous zero locus.
For strictly positive first Chern coordinates, enumerate every nonnegative
ambient curve class with first Chern degree at most `dimension + 1`, including
zero. Each owned curve vector uses compact `marked_nodes` order. Results are
ordered by total coordinate sum, then lexicographically. The index is the gcd
of the ambient first Chern coordinates, rather than their vector or necessarily
the intrinsic Fano index of a zero locus.
If all coordinates vanish, return `(0, [zero_class])`. Mixed zero and positive, or
negative, coordinates require explicit `betas` in
[`small_quantum_multiplication`](@ref).

# Examples
```jldoctest
julia> fano_curve_classes(projective_space(2))
(3, [[0], [1]])

julia> index, betas = fano_curve_classes(product(projective_space(1), projective_space(2)));

julia> index
1

julia> [1, 0] in betas && [0, 1] in betas
true
```
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
target's first Chern class. Return a named tuple with these fields:

- `matrix`: a square matrix of exact [`NovikovPolynomial`](@ref)s. Columns are
  input classes and rows are output coefficients in `basis` order.
- `grading`: a diagonal matrix with entries `codimension - dimension/2`.
- `basis`: ambient Schubert representatives whose restrictions form the chosen basis.
- `metric`: their exact intersection pairing on the target.

All returned arrays, basis words, and polynomial dictionaries are owned by
the caller. [`evaluate`](@ref) specializes Novikov variables; broadcasting
`evaluate.(result.matrix)` sets every variable to one.

For a zero locus, ambient Schubert restrictions are reduced modulo the radical
of the intersection pairing. The ambient subspace need not be preserved by the
full quantum product; nonambient cohomology is excluded. Grading is complex
codimension minus half the target dimension.

Explicit `betas` define a finite truncation: include zero for the classical
term; duplicates are ignored. Otherwise [`fano_curve_classes`](@ref) supplies
the finite list. The defining bundle must satisfy the supported Euler-twist
weight condition, or be a `ConvexEulerTwist` with a geometric convexity
certificate. In the latter case the supplied graded class is used only to
compute the full signed index of the actual bundle.

Set `algorithm=:branches` to resum rooted branches by power-series recursion;
the default `:trees` uses explicit graph enumeration. `progress` may be a
callback receiving degree/state information from the branch backend.

Set `arithmetic=:certified` to use the branch backend over independent prime
fields, with proved denominator and magnitude bounds and deterministic rational
recovery. `workers` selects a bounded Distributed process pool, `checkpoint`
selects an atomic TOML checkpoint directory, and `torus_seeds` selects two
distinct nonsingular exact specializations. The result additionally contains
`certification` records. See [`certified_pair_invariants`](@ref) for checkpoint
semantics. `worker_heap_size` is a per-worker GC hint, not a hard memory limit.

# Examples
```jldoctest
julia> result = small_quantum_multiplication(projective_space(1));

julia> evaluate.(result.matrix) == [0 2; 2 0]
true

julia> evaluate(result.matrix[1, 2], (3,))
6//1
```
"""
function small_quantum_multiplication(X::PartialFlagVariety; bundle=nothing, betas=nothing,
  algorithm=:trees, progress=nothing, arithmetic=:rational, workers=0,
  worker_heap_size=nothing, checkpoint=nothing, torus_seeds=(0, 3), max_primes=256)
  _small_quantum_multiplication(_gw_data(X), bundle, betas; algorithm, progress,
    arithmetic, workers, worker_heap_size, checkpoint, torus_seeds, max_primes)
end

"""Validate target and curve inputs before constructing intersection matrices."""
function _small_quantum_multiplication(data::_GWData, bundle, betas;
  algorithm=:trees, progress=nothing, arithmetic=:rational, workers=0,
  worker_heap_size=nothing, checkpoint=nothing, torus_seeds=(0, 3), max_primes=256)
  _check_localization_algorithm(algorithm)
  arithmetic in (:rational, :certified) ||
    throw(ArgumentError("unknown arithmetic: $arithmetic"))
  if arithmetic == :rational
    (workers == 0 && checkpoint === nothing) ||
      throw(ArgumentError("workers and checkpoints require arithmetic=:certified"))
  end
  dim, c, weights = _complete_intersection_data(data, bundle)
  beta_list = if betas === nothing
    last(_fano_curve_classes(data, dim, c))
  else
    unique([_normalize_beta(data, beta) for beta in betas])
  end
  if arithmetic == :certified
    certificates = Any[]
    result = _cert_with_pool(workers; worker_heap_size) do pool
      pair_evaluator =
        (classes, pairs, beta) -> begin
          result = _cert_pair_block(data, classes, pairs, beta, weights, pool;
            checkpoint, torus_seeds, progress, max_primes)
          push!(certificates, result.certificate)
          result.values
        end
      _quantum_multiplication(Val(picard_rank(data.variety)), data, dim, c, weights,
        beta_list; algorithm=:branches, progress, pair_evaluator)
    end
    return merge(result, (; certification=certificates))
  end
  _quantum_multiplication(Val(picard_rank(data.variety)), data, dim, c, weights,
    beta_list; algorithm, progress)
end

"""
    _quantum_multiplication(Val(P), data, dim, c, weights, beta_list)

Assemble multiplication degree by degree. One exact LU factorization raises
an index of each rational pairing matrix; sparse polynomial terms are added
only after this solve. This avoids a matrix inverse and repeated polynomial
copies in a cubic matrix-multiplication loop. `P` specializes the Novikov rank.
"""
function _quantum_multiplication(::Val{P}, data::_GWData, dim::Int, c, weights,
  beta_list; algorithm=:trees, progress=nothing, pair_evaluator=nothing) where {P}
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
        data.codimensions[indices[i]] > 0 && data.codimensions[indices[j]] > 0 &&
          data.codimensions[indices[i]] + data.codimensions[indices[j]] == required_degree
      ]
      isempty(pairs) && continue
      invariants = if pair_evaluator === nothing
        _checked_pair_block(data, indices, pairs, beta, weights; algorithm, progress)
      else
        pair_evaluator(indices, pairs, beta)
      end
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
  grading = Diagonal([(2 * data.codimensions[i] - dim)//2 for i in indices])
  (; matrix, grading, basis, metric)
end

fano_curve_classes(Z::ZeroLocus) =
  fano_curve_classes(ambient_variety(Z); bundle=defining_bundle(Z))

function small_quantum_multiplication(Z::ZeroLocus; betas=nothing,
  algorithm=:trees, progress=nothing, kwargs...)
  small_quantum_multiplication(ambient_variety(Z); bundle=defining_bundle(Z), betas,
    algorithm, progress, kwargs...)
end
