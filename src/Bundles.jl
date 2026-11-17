"""
    ConvexEulerTwist(graded; certificate)

Euler twist by the equivariant K-class `graded` of an actual genus-zero
convex bundle. `graded` is a `CompletelyReducibleBundle` representing that
K-class, and need not itself be convex or admit a regular section.

The caller supplies a nonempty explanation in `certificate` establishing
the actual bundle and its convexity. To interpret the result as zero-locus
invariants, a regular section of the actual bundle is also required.
These geometric assertions are not inferred from the graded pieces.

Unlike the ordinary bundle input, this input permits negative graded
splittings: the full edge index H^0-H^1 is assembled and equal characters
are cancelled before torus specialization. A remaining negative index
multiplicity is rejected. Identically zero fiber weights still require
fiber-equivariant regularization and are unsupported.
"""
struct ConvexEulerTwist
  graded::CompletelyReducibleBundle
  certificate::String

  function ConvexEulerTwist(graded::CompletelyReducibleBundle;
    certificate::AbstractString)
    isempty(strip(certificate)) &&
      throw(ArgumentError("a convexity certificate is required"))
    new(graded, String(certificate))
  end
end

"""First Chern weight of the flag tangent bundle, in ambient fundamental weights."""
function _chern_coordinates(data::_GWData)
  DT = dynkin_type(data.variety)
  c = WeightLatticeElem(DT)
  for i in data.nonparabolic_indices
    c += WeightLatticeElem(data.positive_roots[i])
  end
  Int.(coefficients(c))
end

"""First Chern weight of a bundle, obtained from its determinant representation."""
function _bundle_chern_coordinates(E::CompletelyReducibleBundle)
  λ = p_dominant_weight(only(components(det(E))))
  Int.(coefficients(λ))
end

"""
    _complete_intersection_data(data, bundle)

Validate a target and return `(dimension, chern_coordinates, weights)`.
`nothing` denotes the untwisted ambient variety; an empty weight vector denotes
a zero-rank defining bundle. These have the same mathematical effect.

Ordinary bundle input requires nonnegative weight–coroot pairings, a
sufficient condition rather than a general convexity test. ConvexEulerTwist
instead carries an explicit geometric convexity certificate for the actual
bundle represented by the supplied K-class.
"""
function _complete_intersection_data(data::_GWData{DT,R}, E) where {DT,R}
  E === nothing && return (dimension(data.variety), _chern_coordinates(data), nothing)
  certified = E isa ConvexEulerTwist
  certified && (E = E.graded)
  E isa CompletelyReducibleBundle || throw(
    ArgumentError("the defining bundle must be a CompletelyReducibleBundle")
  )
  variety(E) == data.variety || throw(ArgumentError("bundle and flag variety do not match"))
  dim = dimension(data.variety) - rank(E)
  dim >= 0 ||
    throw(ArgumentError("the defining bundle has rank greater than the ambient dimension"))
  weights = _bundle_weights(DT, E)
  any(iszero, weights) && throw(
    ArgumentError(
      "a defining bundle with a zero torus weight requires fiber-equivariant regularization"
    ),
  )
  if !certified
    for (q, root_index) in enumerate(data.nonparabolic_indices), λ in weights
      degree_on_curve = _root_coroot_pairing(λ, data.positive_coroots[root_index])
      degree_on_curve >= 0 || throw(
        ArgumentError(
          "the Euler-twist weight condition fails: pairing $degree_on_curve on invariant curve $q"
        ),
      )
    end
  end
  c = _chern_coordinates(data) - _bundle_chern_coordinates(E)
  (dim, c, weights)
end

"""Evaluate a root or weight against a coroot using fundamental-weight coordinates."""
_root_coroot_pairing(root::RootSpaceElem, coroot::RootSpaceElem) =
  _root_coroot_pairing(WeightLatticeElem(root), coroot)

function _root_coroot_pairing(weight::WeightLatticeElem, coroot::RootSpaceElem)
  λ = coefficients(weight)
  sum(λ[i] * coefficients(coroot)[i] for i in eachindex(λ))
end

"""
    _bundle_weights(E)

Expand the Levi components of `E` with multiplicity, retaining a concrete
ambient weight type. Weights use the Borel–Weil–Bott convention of
PartialFlagVarieties; their negatives are the actual fiber characters.
"""
function _bundle_weights(E::CompletelyReducibleBundle)
  _bundle_weights(dynkin_type(variety(E)), E)
end

"""Specialize the weight buffer on the runtime ambient Dynkin type."""
function _bundle_weights(::Type{DT}, E::CompletelyReducibleBundle) where {DT<:DynkinType}
  result = WeightLatticeElem{DT,rank(DT)}[]
  sizehint!(result, rank(E))
  for component in components(E)
    append!(result, _component_weights(component))
  end
  result
end

"""
    _component_weights(component)

Lift each Levi weight to the ambient weight lattice. Differences from the
highest Levi weight are first expressed in Levi simple roots, then inserted
at `levi_permutation` positions in the ambient root system. Freudenthal
multiplicities are preserved.
"""
function _component_weights(component::IrrepLevi)
  mdt = marked_dynkin_type(component)
  λ = p_dominant_weight(component)
  is_borel(mdt) && return [λ]
  ss = semisimple_part(component)
  iszero(ss) && return [λ]

  LT = levi_type(mdt)
  Cinv = cartan_matrix_inverse(LT)
  order = levi_permutation(mdt)
  DT = dynkin_type(mdt)
  weights = typeof(λ)[]
  for (ν, multiplicity) in freudenthal_formula(ss)
    d = Cinv * (coefficients(ss) - ν)
    all(isinteger, d) || error("nonintegral Levi root coordinates")
    root_coefficients = zeros(Int, rank(DT))
    for i in eachindex(order)
      root_coefficients[order[i]] = Int(d[i])
    end
    ambient = λ - WeightLatticeElem(RootSpaceElem(DT, root_coefficients))
    append!(weights, Iterators.repeated(ambient, Int(multiplicity)))
  end
  weights
end
