"""
    _BilleyPlan

Finite-state recurrence for one Schubert restriction. States are prefixes of
reduced expressions of the requested Weyl element, ordered by Coxeter length.
`steps[i, s]` is the state reached by appending generator `s`, or zero if that
product cannot occur in a reduced expression of the target. The first state
is the identity and the last state is the target.
"""
struct _BilleyPlan
  steps::Matrix{Int}
end

"""
    _GWData{DT,R}

Cached, torus-independent geometry of a flag variety. Fixed points use minimal
right-coset representatives. `steps[color, root]` is the opposite endpoint of
the invariant curve indexed by a nonparabolic positive root; `curve_classes`
records its marked-node coroot coordinates.

Root actions, reduced-word roots, neighbors, inverse edge lookups, and
root–coroot pairings are computed once. The only mutable mathematical cache
is `billey_plans`, protected by its own lock. Treat all other stored arrays
and Weyl words as read-only; public Schubert constructors return owned words.
"""
struct _GWData{DT,R}
  variety::PartialFlagVariety
  weyl_group::WeylGroup{DT,R}
  representatives::Vector{WeylGroupElem{DT,R}}
  representative_index::Dict{SVector{R,Int},Int}
  stabilizer_weight::WeightLatticeElem{DT,R}
  positive_roots::Vector{RootSpaceElem{DT,R}}
  positive_coroots::Vector{RootSpaceElem{DT,R}}
  nonparabolic_indices::Vector{Int}
  curve_classes::Vector{Vector{Int}}
  steps::Matrix{Int}
  generators::Vector{WeylGroupElem{DT,R}}
  inverse_representatives::Vector{WeylGroupElem{DT,R}}
  codimensions::Vector{Int}
  root_coordinates::Matrix{SVector{R,Int}}
  word_roots::Vector{Vector{SVector{R,Int}}}
  root_pairings::Matrix{Int}
  neighbors::Vector{Vector{Int}}
  edge_roots::Vector{Dict{Int,Int}}
  billey_plans::Dict{Int,_BilleyPlan}
  billey_lock::ReentrantLock
end

"""Shared root data keyed by marked Dynkin type; heterogeneous values end at typed barriers."""
const _GW_DATA_CACHE = Dict{MarkedDynkinType,Any}()
"""Protect root-data construction and cache clearing; arithmetic uses computation-local state."""
const _GW_DATA_LOCK = ReentrantLock()

"""Return shared flag data, constructing it once under the structural-cache lock."""
function _gw_data(X::PartialFlagVariety)
  mdt = marked_dynkin_type(X)
  lock(_GW_DATA_LOCK) do
    get!(_GW_DATA_CACHE, mdt) do
      _build_gw_data(dynkin_type(X), X)
    end
  end
end

"""Build root-specialized GKM data; the Dynkin type is a specialization barrier."""
function _build_gw_data(::Type{DT}, X::PartialFlagVariety) where {DT<:DynkinType}
  R = rank(DT)
  W = weyl_group(DT)
  reps = right_coset_reps(W, collect(unmarked_nodes(X)))
  inverses = inv.(reps)
  generators = [gen(W, i) for i in 1:R]
  stabilizer = sum(
    (fundamental_weight(DT, i) for i in marked_nodes(X));
    init=WeightLatticeElem(DT),
  )
  rep_index = Dict{SVector{R,Int},Int}(
    coefficients(stabilizer * w_inverse) => i for (i, w_inverse) in enumerate(inverses)
  )

  RS = RootSystem(DT)
  roots = positive_roots(RS)
  coroots = positive_coroots(RS)
  marked = marked_nodes(X)
  nonparabolic = Int[
    i for i in eachindex(roots) if any(!iszero(coefficients(roots[i])[j]) for j in marked)
  ]
  curves = [Int[coefficients(coroots[i])[j] for j in marked] for i in nonparabolic]

  steps = Matrix{Int}(undef, length(reps), length(nonparabolic))
  for (i, w_inverse) in enumerate(inverses), (q, root_index) in enumerate(nonparabolic)
    endpoint = reflect(stabilizer, roots[root_index]) * w_inverse
    steps[i, q] = rep_index[coefficients(endpoint)]
  end

  # The weight stabilizer identifies cosets without enumerating the full Weyl
  # group. A per-vertex inverse lookup also avoids rebuilding tuple-keyed
  # dictionaries for every free tree and every curve class.
  edge_roots = [Dict(steps[i, q] => q for q in axes(steps, 2)) for i in eachindex(reps)]
  neighbors = [sort!(collect(keys(lookup))) for lookup in edge_roots]
  coordinates = [coefficients(root * w_inverse) for w_inverse in inverses, root in roots]
  word_roots = map(reps) do w
    prefix = one(W)
    result = SVector{R,Int}[]
    for s in word(w)
      push!(result, coefficients(simple_root(RS, Int(s)) * inv(prefix)))
      prefix *= generators[Int(s)]
    end
    result
  end
  pairings = [
    _root_coroot_pairing(roots[i], coroots[j]) for i in nonparabolic, j in nonparabolic
  ]
  _GWData{DT,R}(
    X, W, reps, rep_index, stabilizer, roots, coroots, nonparabolic, curves, steps,
    generators, inverses, length.(reps), coordinates, word_roots, pairings,
    neighbors, edge_roots, Dict{Int,_BilleyPlan}(), ReentrantLock(),
  )
end

"""
    _billey_plan(data, class_index)

Return the cached reduced-subword automaton for a Schubert class. Repeatedly
deleting right descents enumerates precisely the prefixes of reduced words
of the target. Other Weyl elements cannot contribute to its restriction and
are excluded before any rational arithmetic is performed.
"""
function _billey_plan(data::_GWData, class_index::Int)
  lock(data.billey_lock) do
    get!(data.billey_plans, class_index) do
      target = data.representatives[class_index]
      states = [target]
      seen = Set(states)
      for w in states, descendant in bruhat_descendants(w)
        descendant in seen && continue
        push!(seen, descendant)
        push!(states, descendant)
      end
      sort!(states; by=length)
      index = Dict(w => i for (i, w) in enumerate(states))
      steps = zeros(Int, length(states), length(data.generators))
      for (i, w) in enumerate(states), (s, reflection) in enumerate(data.generators)
        successor = w * reflection
        length(successor) == length(w) + 1 || continue
        steps[i, s] = get(index, successor, 0)
      end
      _BilleyPlan(steps)
    end
  end
end

"""Copy a Weyl word at the public boundary so callers cannot mutate cached classes."""
_owned_schubert(data::_GWData, w::WeylGroupElem) = data.weyl_group(word(w); normalize=false)

"""
    schubert_classes(X::PartialFlagVariety)

Return the Schubert basis of `X`, represented by the minimal-length elements
of the Weyl-group quotient.  The list is ordered by length and starts with the
identity class.
"""
function schubert_classes(X::PartialFlagVariety)
  data = _gw_data(X)
  [_owned_schubert(data, w) for w in data.representatives]
end

"""Find the minimal representative using the orbit of the parabolic stabilizer weight."""
function _project_to_flag(data::_GWData{DT,R}, w::WeylGroupElem{DT,R}) where {DT,R}
  key = coefficients(data.stabilizer_weight * inv(w))
  data.representatives[data.representative_index[key]]
end

"""
    schubert_class(X::PartialFlagVariety, word)

Return the Schubert class indexed by `word`, projected to the minimal
representative of the corresponding Weyl-group coset.
"""
function schubert_class(X::PartialFlagVariety, word::AbstractVector{<:Integer})
  data = _gw_data(X)
  all(s -> 1 <= s <= length(data.generators), word) ||
    throw(
      ArgumentError("simple-reflection indices must lie in 1:$(length(data.generators))")
    )
  _owned_schubert(data, _project_to_flag(data, data.weyl_group(Vector{Int}(word))))
end

schubert_class(X::PartialFlagVariety, word::Tuple{Vararg{Integer}}) =
  schubert_class(X, collect(word))

"""
    point_class(X::PartialFlagVariety)

Return the top Schubert class of `X`.
"""
function point_class(X::PartialFlagVariety)
  data = _gw_data(X)
  _owned_schubert(data, last(data.representatives))
end

"""Validate and copy a compact or ambient curve class into marked-node order."""
function _normalize_beta(data::_GWData, beta)
  beta isa Union{AbstractVector,Tuple} ||
    throw(ArgumentError("a curve class must be a vector or tuple of integers"))
  all(value -> value isa Integer, beta) ||
    throw(ArgumentError("curve-class coordinates must be integers"))
  all(value -> 0 <= value <= typemax(Int), beta) ||
    throw(ArgumentError("curve-class coordinates must be nonnegative machine integers"))
  values = Int.(collect(beta))
  marked = marked_nodes(data.variety)
  R = rank(dynkin_type(data.variety))
  if length(values) == length(marked)
    compact = values
  elseif length(values) == R
    unmarked_nonzero = [
      (i, values[i]) for i in unmarked_nodes(data.variety) if !iszero(values[i])
    ]
    isempty(unmarked_nonzero) || throw(
      ArgumentError(
        "curve class may be nonzero only at marked nodes; got $unmarked_nonzero"
      ),
    )
    compact = values[collect(marked)]
  else
    throw(
      ArgumentError(
        "curve class needs $(length(marked)) marked-node coordinates or $R ambient coordinates"
      ),
    )
  end
  all(>=(0), compact) || throw(ArgumentError("curve-class coordinates must be nonnegative"))
  compact
end
