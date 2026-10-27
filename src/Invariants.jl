"""
    _curve_degree(data, chern_coordinates, beta)

Pair an ambient first-Chern weight with compact, validated marked-node curve
coordinates. Unmarked coordinates are absent from `beta`.
"""
function _curve_degree(data::_GWData, c, beta)
  degree = 0
  for (node, b) in zip(marked_nodes(data.variety), beta)
    degree = Base.Checked.checked_add(degree, Base.Checked.checked_mul(c[node], b))
  end
  degree
end

"""Maximum fixed-tree vertex count, checking machine-integer overflow before allocation."""
_vertex_bound(beta) = foldl(Base.Checked.checked_add, beta; init=1)

"""Expected complex codimension, including the degree-zero classical pairing convention."""
function _expected_degree(data::_GWData, beta, n::Int, dim::Int, c)
  if all(iszero, beta)
    return Base.Checked.checked_add(dim, max(n - 3, 0))
  end
  Base.Checked.checked_add(
    Base.Checked.checked_add(dim, _curve_degree(data, c, beta)), n - 3
  )
end

"""
    expected_gromov_witten_degree(X, beta, n; bundle=nothing)
    expected_gromov_witten_degree(Z::ZeroLocus, beta, n)

Return the required sum of Schubert codimensions for a genus-zero primary
invariant in curve class `beta` with `n` insertions.  With `bundle=E`, compute
the requirement for the smooth zero locus of the homogeneous bundle `E`,
subject to the supported weight condition. For positive degree the result is
`dim + c1⋅beta + n - 3`.
At degree zero, fewer than three marks use the classical intersection
convention; otherwise the expected degree is `dim + n - 3`.
"""
function expected_gromov_witten_degree(
  X::PartialFlagVariety, beta, n::Integer; bundle=nothing
)
  n >= 0 || throw(ArgumentError("the number of insertions must be nonnegative"))
  data = _gw_data(X)
  β = _normalize_beta(data, beta)
  dim, c, _ = _complete_intersection_data(data, bundle)
  _expected_degree(data, β, Int(n), dim, c)
end

"""
    gromov_witten_invariant(X, classes, beta; bundle=nothing)
    gromov_witten_invariant(Z::ZeroLocus, classes, beta)

Compute an exact genus-zero primary Gromov--Witten invariant of the flag
variety `X`.  `classes` are entries returned by [`schubert_classes`](@ref), and
`beta` is written either in marked-node order or as an ambient zero-padded
vector. Passing a homogeneous `bundle` satisfying the supported weight
condition computes the Euler-twisted invariant of its smooth zero locus.
Every defining weight must pair nonnegatively with every nonparabolic
positive coroot, and zero weights are unsupported. Smoothness and regularity
of the defining section are assumptions. The result is a `Rational{BigInt}`.
Dimension-mismatched invariants and degree-zero primary invariants with
more than three marks vanish. With at most three marks, degree zero uses
the classical intersection integral, including the two-point pairing.
"""
function gromov_witten_invariant(
  X::PartialFlagVariety, classes::AbstractVector, beta; bundle=nothing
)
  _gromov_witten_invariant(_gw_data(X), classes, beta, bundle)
end

"""Validate insertions and the dimension axiom before allocating any torus evaluation."""
function _gromov_witten_invariant(data::_GWData, classes, beta, bundle)
  β = _normalize_beta(data, beta)
  class_indices = Int[_class_index(data, class) for class in classes]
  dim, c, weights = _complete_intersection_data(data, bundle)
  expected = _expected_degree(data, β, length(classes), dim, c)
  sum((data.codimensions[i] for i in class_indices); init=0) == expected ||
    return zero(_BigRational)
  # For more than three marks, primary insertions have no degree on M̄₀,ₙ.
  all(iszero, β) && length(classes) > 3 && return zero(_BigRational)
  _checked_localization(data, class_indices, β, weights)
end

function expected_gromov_witten_degree(Z::ZeroLocus, beta, n::Integer)
  expected_gromov_witten_degree(
    ambient_variety(Z), beta, n; bundle=defining_bundle(Z)
  )
end

function gromov_witten_invariant(Z::ZeroLocus, classes::AbstractVector, beta)
  gromov_witten_invariant(
    ambient_variety(Z), classes, beta; bundle=defining_bundle(Z)
  )
end
