"""
    LogMonodromyBlock

One eigenvalue group of a local monodromy, represented by its rational exponent
in `[0,1)` and the Jordan block sizes of its normalized logarithm. The actual
monodromy eigenvalue is `exp(2πim * exponent)`. `ranks` records the ranks of
successive powers of the nilpotent logarithm, starting with the identity.
"""
struct LogMonodromyBlock
  exponent::QQFieldElem
  sizes::Vector{Int}
  ranks::Vector{Int}
  resonance_width::Int
end

"""
    LocalMonodromyData

Exact complex local conjugacy data. `locus` is `:zero`, `:infinity`, or a monic
irreducible polynomial over `Nemo.QQ`; in the last case the data applies at each
of its `degree` roots. `exponents` are the indicial roots before reduction modulo
integers. `fixed_dimension` is the dimension of the monodromy invariant space,
and `ramification` is its codimension at **one** root. `apparent` means that the
local monodromy is the identity, including at an ordinary point.
"""
struct LocalMonodromyData
  locus::Union{Symbol,QQPolyRingElem}
  degree::Int
  exponents::Vector{QQFieldElem}
  blocks::Vector{LogMonodromyBlock}
  fixed_dimension::Int
  ramification::Int
  apparent::Bool
end

"""Convert exact coefficient rows to polynomials, remove their common factor, and validate order."""
function _monodromy_polynomials(coefficients)
  R, _ = polynomial_ring(QQ, :t)
  p = QQPolyRingElem[]
  for row in coefficients
    values = row isa QQPolyRingElem ? collect(Nemo.coefficients(row)) : collect(row)
    all(x -> x isa Integer || x isa Rational || x isa QQFieldElem, values) ||
      throw(ArgumentError("operator coefficients must be exact rational numbers"))
    push!(p, R(QQ.(values)))
  end
  length(p) >= 2 && !iszero(last(p)) ||
    throw(ArgumentError("a nonzero leading coefficient and positive order are required"))
  common = foldl(gcd, p)
  return [divexact(v, common) for v in p]
end

"""Return all roots with multiplicity, rejecting indicial polynomials without rational roots."""
function _rational_indicial_roots(p::QQPolyRingElem, n::Int)
  degree(p) == n || throw(ArgumentError("the point is not regular singular"))
  result = QQFieldElem[]
  for (f, multiplicity) in factor(p)
    degree(f) == 1 || throw(ArgumentError("rational local exponents are required"))
    append!(result, fill(-coeff(f, 0) / coeff(f, 1), multiplicity))
  end
  return sort!(result)
end

"""Compute exact resonant Frobenius kernels over the coefficient field; rows index powers of log."""
function _frobenius_blocks(local_coefficients, exponents; max_resonance::Int)
  n = length(local_coefficients) - 1
  K = parent(local_coefficients[end][1])
  R, z = polynomial_ring(K, :lambda)
  residues = sort!(unique([e - floor(BigInt, e) for e in exponents]))
  result = LogMonodromyBlock[]
  for residue in residues
    group = filter(e -> e - floor(BigInt, e) == residue, exponents)
    lowest = minimum(group)
    width = Int(numerator(maximum(group) - lowest))
    width <= max_resonance || throw(ArgumentError("resonance width exceeds max_resonance"))
    size = (width + 1) * n
    relations = zero_matrix(K, size, size)
    for k in 0:width
      polynomial = R([k < length(row) ? row[k + 1] : zero(K) for row in local_coefficients])
      for step in k:width
        shifted = Nemo.evaluate(polynomial, z + K(lowest + step - k))
        for i in 0:(n - 1), j in i:(n - 1)
          relations[step * n + i + 1, (step - k) * n + j + 1] = coeff(shifted, j - i)
        end
      end
    end
    dimension, solutions = nullspace(relations)
    dimension == length(group) || error("incorrect resonant Frobenius solution dimension")
    ranks = [dimension]
    # In the divided-power basis log(s)^j/j!, differentiation lowers j.
    # Past the last indicial root all further coefficients extend uniquely.
    for power in 1:(n + 1)
      transformed = zero_matrix(K, size, dimension)
      for step in 0:width, j in 0:(n - power - 1), column in 1:dimension
        transformed[step * n + j + 1, column] = solutions[step * n + j + power + 1, column]
      end
      push!(ranks, Nemo.rank(transformed))
    end
    sizes = Int[]
    for block_length in 1:n
      count = ranks[block_length] - 2ranks[block_length + 1] + ranks[block_length + 2]
      count >= 0 || error("invalid nilpotent rank sequence")
      append!(sizes, fill(block_length, count))
    end
    sum(sizes) == length(group) || error("incorrect Jordan block dimensions")
    push!(result, LogMonodromyBlock(residue, sort!(sizes; rev=true), ranks, width))
  end
  return result
end

"""Assemble invariant dimensions, counting only logarithm blocks with eigenvalue one."""
function _local_monodromy_record(locus, degree, exponents, blocks)
  fixed = sum((length(b.sizes) for b in blocks if iszero(b.exponent)); init=0)
  ramification = length(exponents) - fixed
  return LocalMonodromyData(
    locus, degree, exponents, blocks, fixed, ramification, ramification == 0
  )
end

"""Convert Euler-derivative coefficients to ordinary derivative coefficients using Stirling numbers."""
function _ordinary_derivative_coefficients(p)
  n = length(p) - 1
  t = Nemo.gen(parent(p[1]))
  stirling = zeros(BigInt, n + 1, n + 1)
  stirling[1, 1] = 1
  for m in 1:n, k in 1:m
    stirling[m + 1, k + 1] = k * stirling[m, k + 1] + stirling[m, k]
  end
  return [
    sum((stirling[m + 1, k + 1] * p[m + 1] for m in k:n); init=zero(t)) * t^k for k in 0:n
  ]
end

"""Return the derivative divided by its factorial, avoiding machine-integer factorial overflow."""
function _taylor_derivative(p, k::Int)
  for j in 1:k
    p = derivative(p) / j
  end
  return p
end

"""Construct a block record from known exact Jordan sizes, including its nilpotent ranks."""
function _known_logarithm_block(exponent, sizes, width, n)
  ranks = [sum(max(s - power, 0) for s in sizes) for power in 0:(n + 1)]
  return LogMonodromyBlock(exponent, sizes, ranks, width)
end

"""Handle simple leading zeros with a repeated exponent or the exact first Taylor obstruction."""
function _simple_zero_blocks(ordinary, f, exponents)
  n = length(exponents)
  exceptional = copy(exponents)
  for j in 0:(n - 2)
    position = findfirst(==(QQ(j)), exceptional)
    isnothing(position) && return nothing
    deleteat!(exceptional, position)
  end
  sigma = only(exceptional)
  residue = sigma - floor(BigInt, sigma)
  n == 1 && return [_known_logarithm_block(residue, [1], 0, n)]
  if !iszero(residue)
    return [
      _known_logarithm_block(QQ(0), ones(Int, n - 1), max(n - 2, 0), n),
      _known_logarithm_block(residue, [1], 0, n),
    ]
  end
  repeated = 0 <= sigma <= n - 2
  if repeated
    sizes = vcat([2], ones(Int, n - 2))
  elseif sigma == n && n >= 2
    # Eliminate the first Taylor coefficient equation, then clear q_n'(a).
    # The second equation is the only remaining resonant compatibility test.
    leading = rem(derivative(ordinary[end]), f)
    combination = rem(derivative(ordinary[end - 1]) + ordinary[end - 2], f)
    apparent = all(0:(n - 2)) do j
      previous = j == 0 ? zero(f) : ordinary[j]
      iszero(
        rem(
          combination * ordinary[j + 1] +
          leading * (derivative(ordinary[j + 1]) + previous),
          f,
        ),
      )
    end
    sizes = apparent ? ones(Int, n) : vcat([2], ones(Int, n - 2))
  else
    return nothing
  end
  width = Int(numerator(maximum(exponents) - minimum(exponents)))
  return [_known_logarithm_block(QQ(0), sizes, width, n)]
end

"""Compute the indicial polynomial at an algebraic root without constructing a splitting field."""
function _finite_indicial(ordinary, f, multiplicity)
  n = length(ordinary) - 1
  values = QQPolyRingElem[]
  for k in 0:n
    required = multiplicity - n + k
    value = ordinary[k + 1]
    if required > 0 && !iszero(value)
      iszero(rem(value, f^required)) ||
        throw(ArgumentError("the point is not regular singular"))
    end
    push!(values, required < 0 ? zero(f) : rem(_taylor_derivative(value, required), f))
  end
  leading = last(values)
  position = findfirst(j -> !iszero(coeff(leading, j)), 0:degree(leading))
  isnothing(position) && error("zero leading indicial coefficient")
  index = position - 1
  ratios = [coeff(value, index) / coeff(leading, index) for value in values]
  all(value == ratio * leading for (value, ratio) in zip(values, ratios)) ||
    throw(ArgumentError("a rational indicial polynomial is required"))
  R, z = polynomial_ring(QQ, :lambda)
  falling = one(R)
  indicial = zero(R)
  for k in 0:n
    k > 0 && (falling *= z - (k - 1))
    indicial += ratios[k + 1] * falling
  end
  return _rational_indicial_roots(indicial, n)
end

"""Compute the local Frobenius coefficients over a single singular-point number field."""
function _finite_monodromy(p, f, multiplicity; max_resonance::Int)
  ordinary = _ordinary_derivative_coefficients(p)
  exponents = _finite_indicial(ordinary, f, multiplicity)
  if multiplicity == 1
    blocks = _simple_zero_blocks(ordinary, f, exponents)
    isnothing(blocks) || return _local_monodromy_record(f, degree(f), exponents, blocks)
  end
  n = length(p) - 1
  residues = unique([e - floor(BigInt, e) for e in exponents])
  width = maximum(residues) do residue
    group = filter(e -> e - floor(BigInt, e) == residue, exponents)
    Int(numerator(maximum(group) - minimum(group)))
  end
  width <= max_resonance || throw(ArgumentError("resonance width exceeds max_resonance"))
  K, a = degree(f) == 1 ? (QQ, -coeff(f, 0) / coeff(f, 1)) : number_field(f, :a)
  R, z = polynomial_ring(K, :lambda)
  falling = [one(R)]
  for k in 0:(n - 1)
    push!(falling, falling[end] * (z - k))
  end
  local_polynomials = [zero(R) for _ in 0:width]
  for j in 0:width, k in 0:n
    derivative_order = j + multiplicity - n + k
    derivative_order < 0 && continue
    value = rem(_taylor_derivative(ordinary[k + 1], derivative_order), f)
    scalar = Nemo.evaluate(polynomial_ring(K, :t)[1](collect(Nemo.coefficients(value))), a)
    local_polynomials[j + 1] += scalar * falling[k + 1]
  end
  coefficients = [[coeff(v, m) for v in local_polynomials] for m in 0:n]
  blocks = _frobenius_blocks(coefficients, exponents; max_resonance)
  return _local_monodromy_record(f, degree(f), exponents, blocks)
end

"""
    local_monodromy(coefficients, point=:zero; max_resonance=100)

Compute exact local logarithmic Jordan data for `L = sum(p_m(t) D^m)` with
`D=t*d/dt`. Input row `m+1` lists the coefficients of `p_m` in ascending powers
of `t`; vectors of `Nemo.QQ` polynomials are also accepted. `point` is `:zero`,
`:infinity`, or an irreducible rational polynomial describing a finite locus.
The operator must be regular singular at the requested point and have rational
local exponents. Unsupported inputs throw `ArgumentError`; exponent coincidence
alone is never used to infer logarithms. `max_resonance` bounds the finite
Frobenius system and throws if exceeded.

The result gives complex local conjugacy classes, not connection matrices in a
shared global or integral basis. The method uses the resonant Frobenius kernel
in divided logarithmic powers, equivalent to the prepared-residue method of
Coates--Galkin--Kasprzyk--Strangeway, arXiv:1406.4891, section 9.

```jldoctest
julia> local_monodromy([[0], [0], [1]]).blocks[1].sizes
1-element Vector{Int64}:
 2

julia> local_monodromy([[0], [-1], [1]]).blocks[1].sizes
2-element Vector{Int64}:
 1
 1
```
"""
function local_monodromy(coefficients, point=:zero; max_resonance::Int=100)
  p = _monodromy_polynomials(coefficients)
  n = length(p) - 1
  if point isa QQPolyRingElem
    degree(point) >= 1 && is_irreducible(point) ||
      throw(ArgumentError("the finite locus must be nonconstant and irreducible"))
    f = parent(p[1])(collect(Nemo.coefficients(point / leading_coefficient(point))))
    f == Nemo.gen(parent(f)) && return local_monodromy(p, :zero; max_resonance)
    multiplicity = 0
    leading = last(p)
    while iszero(rem(leading, f))
      multiplicity += 1
      leading = divexact(leading, f)
    end
    # At an ordinary point the monodromy is the identity.
    if multiplicity == 0
      exponents = QQ.(0:(n - 1))
      blocks = [_known_logarithm_block(QQ(0), ones(Int, n), n - 1, n)]
      return _local_monodromy_record(f, degree(f), exponents, blocks)
    end
    return _finite_monodromy(p, f, multiplicity; max_resonance)
  end
  point in (:zero, :infinity) || throw(ArgumentError("unknown local point"))
  if point == :infinity
    d = maximum(degree, p)
    R = parent(p[1])
    p = [R([(-1)^(i - 1) * coeff(v, d - k) for k in 0:d]) for (i, v) in enumerate(p)]
  end
  R, _ = polynomial_ring(QQ, :lambda)
  exponents = _rational_indicial_roots(R([coeff(v, 0) for v in p]), n)
  rows = [collect(Nemo.coefficients(v)) for v in p]
  # Retain a zero coefficient for identically zero polynomials.
  rows = [isempty(row) ? [QQ(0)] : row for row in rows]
  blocks = _frobenius_blocks(rows, exponents; max_resonance)
  return _local_monodromy_record(point, 1, exponents, blocks)
end

"""
    log_monodromy_matrix(data::LocalMonodromyData)

Return a rational Jordan matrix `B` such that the local monodromy is conjugate
over the complex numbers to `exp(2πim * B)`. Eigenvalue representatives lie in
`[0,1)`. This is the paper's normalized log-monodromy convention.
"""
function log_monodromy_matrix(data::LocalMonodromyData)
  result = zero_matrix(QQ, length(data.exponents), length(data.exponents))
  offset = 0
  for block in data.blocks, size in block.sizes
    for j in 1:size
      result[offset + j, offset + j] = block.exponent
      j < size && (result[offset + j, offset + j + 1] = 1)
    end
    offset += size
  end
  return result
end

"""
    ramification_data(coefficients; max_resonance=100)

Compute local data at zero, every nonzero root of the leading Euler coefficient,
and infinity. Return `(order, local_data, ramification, defect, defect_zero)`.
Ramification sums `rank(T-I)`, weighting each irreducible finite locus by its
degree, and `defect = ramification - 2order`. Apparent singularities contribute
zero. Nontrivial eigenvalues contribute their full generalized eigenspace;
using only the rank of the nilpotent logarithm would give an incorrect answer.

A zero defect is called extremal in arXiv:1406.4891. Its interpretation as the
minimum possible ramification requires a nontrivial irreducible local system;
this function does not assert or test irreducibility or minimality of an input
annihilator for a particular period.
"""
function ramification_data(coefficients; max_resonance::Int=100)
  p = _monodromy_polynomials(coefficients)
  local_data = [local_monodromy(p, :zero; max_resonance)]
  factors = sort!(collect(factor(last(p))); by=x -> (degree(x[1]), string(x[1])))
  for (f, multiplicity) in factors
    iszero(coeff(f, 0)) && continue
    f = f / leading_coefficient(f)
    push!(local_data, _finite_monodromy(p, f, multiplicity; max_resonance))
  end
  push!(local_data, local_monodromy(p, :infinity; max_resonance))
  ramification = sum(point.degree * point.ramification for point in local_data)
  order = length(p) - 1
  defect = ramification - 2order
  return (; order, local_data, ramification, defect, defect_zero=iszero(defect))
end

"""
    ramification_defect(coefficients; kwargs...)

Return `ramification_data(coefficients; kwargs...).defect`. Pass an already
computed ramification result instead to avoid recomputing the local data.
"""
ramification_defect(coefficients; kwargs...) =
  ramification_data(coefficients; kwargs...).defect

"""Read the defect of a previously computed ramification result without recomputation."""
ramification_defect(data::NamedTuple) = data.defect
