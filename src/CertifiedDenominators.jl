"""Positive rational content of a nonzero rational linear form."""
function _cert_content(coordinates)
  scale = foldl(lcm, (denominator(x) for x in coordinates); init=big(1))
  content = foldl(gcd, (BigInt(numerator(x * scale)) for x in coordinates);
    init=big(0))
  iszero(content) && throw(ArgumentError("zero formal character"))
  abs(content)//scale
end

"""LCM of the reduced positive denominators of rational factors."""
_cert_denominator_lcm(values) =
  foldl(lcm, (denominator(x) for x in values); init=big(1))

"""Ordinary prime valuation of a nonzero integer or rational."""
function _cert_valuation(value::Integer, prime)
  iszero(value) && throw(ArgumentError("valuation of zero"))
  count = 0
  while iszero(rem(value, prime))
    count += 1
    value = div(value, prime)
  end
  count
end
_cert_valuation(value::Rational, prime) =
  _cert_valuation(numerator(value), prime) - _cert_valuation(denominator(value), prime)

"""A lower bound for the sum of edge valuations at exactly the requested degree."""
function _cert_edge_valuation_minimum(beta, costs, values)
  strides, size = _branch_strides(beta)
  unreachable = typemax(Int) ÷ 4
  minimums = fill(unreachable, size)
  minimums[1] = 0
  for (index, coordinate) in enumerate(CartesianIndices(Tuple(beta .+ 1)))
    index == 1 && continue
    degree = collect(Tuple(coordinate) .- 1)
    for a in eachindex(costs)
      all(costs[a] .<= degree) || continue
      previous = _branch_degree_index(degree - costs[a], strides)
      minimums[previous] == unreachable && continue
      minimums[index] = min(minimums[index], minimums[previous] + values[a])
    end
  end
  minimums[end] == unreachable && error("requested degree is unreachable")
  minimums
end

"""Gauss content of the normalized edge factor; Weyl conjugation preserves it."""
function _cert_edge_content(evaluation, roots, root, cover)
  data = evaluation.data
  alpha = roots[root]
  a = _cert_content(alpha)
  result = big(cover)^(2 * cover - 1)//(factorial(big(cover))^2) / a^(2 * cover)
  for q in eachindex(roots)
    q == root && continue
    gamma = roots[q]
    m = cover * data.root_pairings[q, root]
    if m <= -2
      for j in 1:(-m - 1)
        result *= _cert_content(gamma + (big(j)//cover) * alpha)
      end
    elseif m >= 0
      for k in 0:m
        result /= _cert_content(gamma - (big(k)//cover) * alpha)
      end
    end
  end
  terms = _build_edge_index_terms(data, evaluation.weight_coordinates,
    evaluation.splittings, root, cover)
  for (j, shift, multiplicity) in terms
    character = -evaluation.weight_coordinates[:, j] + (big(shift)//cover) * alpha
    result *= _cert_content(character)^multiplicity
  end
  result
end

"""
Bound the denominator of the constant equivariant integral in flag degree beta.
Insertions must be integral Schubert classes. The bound is derived from the original graph factors, before torus specialization.
The bound controls the invariant, not the denominators at a numerical torus.
"""
function _cert_denominator_certificate(evaluation, beta; nmarks=2)
  data = evaluation.data
  n = sum(beta)
  n > 0 || throw(ArgumentError("positive curve degree required"))
  roots = [
    BigInt.(
      coefficients(
        WeightLatticeElem(data.positive_roots[index])),
    )
    for index in data.nonparabolic_indices
  ]
  root_contents = _cert_content.(roots)
  fiber_contents = [
    _cert_content(evaluation.weight_coordinates[:, j])
    for j in axes(evaluation.weight_coordinates, 2)
  ]
  euler_content = prod(root_contents) / prod(fiber_contents)
  root_indices, covers, costs = Int[], Int[], Vector{Int}[]
  edge_contents, inverse_contents, leaf_contents = _BigRational[],
  _BigRational[],
  _BigRational[]
  for (root, curve) in enumerate(data.curve_classes)
    maximum = minimum(beta[i] ÷ curve[i] for i in eachindex(beta) if curve[i] > 0)
    for cover in 1:maximum
      push!(root_indices, root)
      push!(covers, cover)
      push!(costs, cover .* curve)
      push!(edge_contents, _cert_edge_content(evaluation, roots, root, cover))
      push!(inverse_contents, cover / root_contents[root])
      push!(leaf_contents, root_contents[root] / cover)
    end
  end
  biv = big(1)
  biv_minimum_cover_sum = Dict{BigInt,Int}()
  for a in eachindex(covers), b in a:length(covers)
    all(costs[a] + costs[b] .<= beta) || continue
    linear = covers[a] * roots[root_indices[b]] + covers[b] * roots[root_indices[a]]
    value = big(covers[a]) * covers[b] / _cert_content(linear)
    biv = lcm(biv, denominator(value))
    den = denominator(value)
    if den > 1
      biv_minimum_cover_sum[den] = min(get(biv_minimum_cover_sum, den, typemax(Int)),
        covers[a] + covers[b])
    end
  end
  factors = Dict(
    "edge" => _cert_denominator_lcm(edge_contents),
    "vertex_euler" => denominator(euler_content),
    "inverse_weight" => _cert_denominator_lcm(inverse_contents),
    "leaf" => _cert_denominator_lcm(leaf_contents), "bivalent" => biv,
  )
  coarse_bound =
    factorial(big(n + 1)) * factors["edge"]^n *
    factors["vertex_euler"]^(n - 1) * factors["inverse_weight"]^(4n + nmarks) *
    factors["leaf"]^(n + 1) * factors["bivalent"]^(n - 1)
  # The sharper strategy applies when the following content checks establish that
  # vertex Euler contents are units and all inverse flag contents are integral.
  if euler_content != 1 || factors["inverse_weight"] != 1
    return Dict("flag_degree" => collect(beta), "marks" => nmarks,
      "edge_bound" => n, "denominator_bound" => string(coarse_bound),
      "denominator_bound_bits" => ndigits(coarse_bound; base=2),
      "strategy" => "coarse graph-factor bound",
      "factor_lcms" => Dict(k => string(v) for (k, v) in factors),
      "basis_hypothesis" => "integral Schubert insertions")
  end
  prime_source = factors["edge"] * factors["leaf"] * biv * factorial(big(n + 1))
  primes = sort([BigInt(p) for (p, _) in Nemo.factor(Nemo.ZZ(prime_source))])
  prime_records = Any[]
  bound = big(1)
  for prime in primes
    edge_values = [
      _cert_valuation(value, prime) - 2 * _cert_valuation(denominator(leaf), prime)
      for (value, leaf) in zip(edge_contents, leaf_contents)
    ]
    minimums = _cert_edge_valuation_minimum(beta, costs, edge_values)
    ratio = maximum(
      (_cert_valuation(den, prime)//count
       for (den, count) in biv_minimum_cover_sum); init=big(0)//1)
    biv_exponent = Int(fld(numerator(2 * n * ratio), denominator(2 * n * ratio)))
    automorphism_exponent = _cert_valuation(factorial(big(n + 1)), prime)
    exponent = max(0, -minimums[end] + biv_exponent + automorphism_exponent)
    bound *= prime^exponent
    push!(
      prime_records,
      Dict("prime" => string(prime), "exponent" => exponent,
        "edge_valuation_minimum" => minimums[end], "all_degree_edge_minima" => minimums,
        "bivalent_valuation_per_cover_sum" => string(ratio),
        "bivalent_exponent" => biv_exponent,
        "automorphism_exponent" => automorphism_exponent),
    )
  end
  @assert iszero(rem(coarse_bound, bound))
  Dict("flag_degree" => beta, "marks" => nmarks,
    "edge_bound" => n, "denominator_bound" => string(bound),
    "denominator_bound_bits" => ndigits(bound; base=2),
    "coarse_denominator_bound_bits" => ndigits(coarse_bound; base=2),
    "prime_bounds" => prime_records,
    "factor_lcms" => Dict(k => string(v) for (k, v) in factors),
    "edge_contents" => [
      Dict("root" => r, "cover" => c, "content" => string(v))
      for (r, c, v) in zip(root_indices, covers, edge_contents)
    ],
    "basis_hypothesis" => "integral Schubert insertions in the ambient flag variety",
    "strategy" => "sharp Gauss valuation bound")
end
