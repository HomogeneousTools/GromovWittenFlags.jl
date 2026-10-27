import GromovWittenFlags as GW
import PartialFlagVarieties as PF
import Semisimple as SS
import GKMtools as GK
import Oscar as OS

"""Common exact scalar field used only when comparing native output matrices."""
const QQ = Rational{BigInt}

"""Construct the same fourfold independently in PartialFlagVarieties and GKMtools."""
function fourfold_models(name::String)
  projective = Dict(
    "cubic" => (5, [3]), "quartic" => (5, [4]),
    "quadric-cubic" => (6, [2, 3]), "three-quadrics" => (7, [2, 2, 2]),
  )
  if haskey(projective, name)
    n, degrees = projective[name]
    X = PF.projective_space(n)
    E = reduce(+, [PF.line_bundle(X, d) for d in degrees])
    return PF.zero_locus(E), GK.vector_bundle_O(n, degrees)
  elseif name == "K3-29"
    X = PF.product(PF.projective_space(2), PF.projective_space(4))
    E = PF.line_bundle(X, [1, 2]) + PF.line_bundle(X, [1, 1])
    return PF.zero_locus(E), GK.vector_bundle_O([2, 4], [[1, 2], [1, 1]])
  end
  k, n = if name == "c5"
    (3, 7)
  elseif name == "b8"
    (2, 7)
  elseif name == "Gushel-Mukai"
    (2, 5)
  else
    throw(ArgumentError("unknown fourfold: $name"))
  end
  X = PF.Gr(k, n)
  S, Q = PF.universal_subbundle(X), PF.universal_quotient_bundle(X)
  Sg, Qg = GK.tautological_and_univ_bd(GK.GKMGraph, k, n)
  Lg = OS.det(OS.dual(Sg))
  E, V = if name == "c5"
    (PF.exterior_power(PF.dual(S), 2) + PF.exterior_power(Q, 3) + PF.line_bundle(X, 1),
      GK.wedge_product(OS.dual(Sg), 2) + GK.wedge_product(Qg, 3) + Lg)
  elseif name == "b8"
    (3 * PF.line_bundle(X, 1) + PF.symmetric_power(PF.dual(S), 2),
      Lg + Lg + Lg + GK.sym_product(OS.dual(Sg), 2))
  else
    (PF.line_bundle(X, 1) + PF.line_bundle(X, 2), Lg + Lg^2)
  end
  PF.zero_locus(E), V
end

"""Convert exact Oscar scalars; reject any unevaluated equivariant parameter."""
as_rational(x::Integer) = QQ(x)
as_rational(x::OS.ZZRingElem) = BigInt(x) // big(1)
as_rational(x::OS.QQFieldElem) = BigInt(numerator(x)) // BigInt(denominator(x))
function as_rational(x)
  n, d = numerator(x), denominator(x)
  OS.is_constant(n) && OS.is_constant(d) || error("expected a nonequivariant scalar")
  as_rational(OS.constant_coefficient(n)) / as_rational(OS.constant_coefficient(d))
end

"""Complete symmetric polynomials h₀,…,hₙ, updated one Chern root at a time."""
function complete_symmetric(roots, n::Int)
  h = [one(first(roots)); fill(zero(first(roots)), n)]
  for root in roots, degree in 1:n
    h[degree + 1] += root * h[degree]
  end
  h
end

"""Jacobi–Trudi formula for an ordinary Schubert class lifted equivariantly via S* roots."""
function schur_value(roots, partition)
  isempty(partition) && return one(first(roots))
  n = length(partition)
  h = complete_symmetric(roots, first(partition) + n - 1)
  entries = [
    let degree = partition[i] - i + j
      degree < 0 ? zero(first(roots)) : h[degree + 1]
    end for i in 1:n, j in 1:n
  ]
  OS.det(OS.matrix(parent(first(roots)), entries))
end

"""
Match ordinary Schubert classes on native type-A flag graphs.

GKMtools' fast Grassmannian constructor uses block permutations, avoiding
Oscar Weyl objects entirely. Reduced words act by adjacent position swaps.
For first-block images i₁<…<iₖ, the partition is (iₖ-k,…,i₁-1).
Billey's formula gives sparse Schubert restrictions by default. Unshifted
Schur polynomials in the Chern roots of S* provide an alternative lift of the
same ordinary classes. The classical metric checks the transport in both cases.
"""
function homogeneous_basis(X, G, classes; basis_lift::Symbol=:billey)
  flags = [v.flag for v in GK.vertices_structure(G)]
  first(flags) isa Tuple || error("expected native type-A block permutations")
  n = length(first(flags))
  k = only(PF.marked_nodes(X))
  variables = collect(OS.gens(GK.equivariant_coefficient_ring(G)))
  partitions = map(classes) do w
    permutation = collect(1:n)
    for s in SS.word(w)
      i = Int(s)
      permutation[i], permutation[i + 1] = permutation[i + 1], permutation[i]
    end
    issorted(permutation[1:k]) && issorted(permutation[(k + 1):end]) ||
      error("the word does not give a minimal right-coset representative")
    partition = reverse([permutation[i] - i for i in 1:k])
    filter!(!iszero, partition)
    sum(partition; init=0) == length(w) || error("word transport changed codimension")
    partition
  end
  basis = if basis_lift == :schur
    [
      GK.polynomial_class(G,
        [schur_value([-variables[flag[i]] for i in 1:k], partition) for flag in flags])
      for partition in partitions
    ]
  elseif basis_lift == :billey
    billey_basis(G, classes, flags, variables)
  else
    throw(ArgumentError("the basis lift must be billey or schur"))
  end
  index = Dict(flag => i for (i, flag) in enumerate(flags))
  simple = collect(1:n)
  simple[k], simple[k + 1] = simple[k + 1], simple[k]
  edge = GK.Edge(index[Tuple(1:n)], index[Tuple(simple)])
  basis, [GK.curve_class(G, edge)], [edge]
end

"""Recover an Oscar Weyl element from an image permutation by undoing adjacent inversions."""
function permutation_to_weyl(flag, W, reflections)
  permutation = collect(flag)
  undo = Int[]
  while !issorted(permutation)
    i = findfirst(j -> permutation[j] > permutation[j + 1], 1:(length(permutation) - 1))
    push!(undo, i)
    permutation[i], permutation[i + 1] = permutation[i + 1], permutation[i]
  end
  result = one(W)
  for i in Iterators.reverse(undo)
    result *= reflections[i]
  end
  result
end

"""
Compute matching Billey lifts with the upstream recurrence, independently of GW's implementation.

The pinned GKMtools public Schubert API requires Weyl-backed vertices, whereas
its fast type-A constructors store permutations. These two pinned private
helpers accept explicit Weyl elements, allowing its own recurrence to evaluate
on the native graph without constructing a slower duplicate generalized graph.
Only input preparation uses them; timed quantum calculations use public APIs.
"""
function billey_basis(G, classes, flags, variables)
  R = OS.root_system(:A, length(first(flags)) - 1)
  W = OS.weyl_group(R)
  reflections = OS.reflection.(OS.simple_roots(R))
  fixed = [permutation_to_weyl(flag, W, reflections) for flag in flags]
  map(classes) do class
    u = one(W)
    for s in SS.word(class)
      u *= reflections[Int(s)]
    end
    length(u) == length(class) || error("word transport changed codimension")
    interval = GK._bruhat_lower_interval(u, reflections)
    values = [GK._billey_localization(G, u, v, variables, interval) for v in fixed]
    GK.polynomial_class(G, values)
  end
end

"""First Chern class of one line summand, using its public fiber-weight accessor."""
function summand_chern(V, summand::Int)
  G = GK.baseof(V)
  variables = collect(OS.gens(GK.equivariant_coefficient_ring(G)))
  values = [
    sum(GK.fiber_weight(V, v, summand)[i] * variables[i]
        for i in eachindex(variables)) for v in 1:GK.num_vertices(G)
  ]
  GK.polynomial_class(G, values)
end

"""Match the two projective-factor classes and primitive curves without using vertex order."""
function product_basis(G, V, classes; basis_lift::Symbol=:billey)
  # The independent native model has summands O(1,2) and O(1,1).
  # Their difference is H₂, and twice the second minus the first is H₁.
  A, B = summand_chern(V, 1), summand_chern(V, 2)
  H1, H2 = 2B - A, A - B
  variables = collect(OS.gens(GK.equivariant_coefficient_ring(G)))
  basis = map(classes) do w
    d1 = count(s -> s <= 2, SS.word(w))
    d2 = length(w) - d1
    if basis_lift == :schur
      H1^d1 * H2^d2
    elseif basis_lift == :billey
      # For Pⁿ the degree-d opposite Schubert class restricts at j to
      # ∏_{i=1}^d(tᵢ-tⱼ). Product restrictions multiply factorwise.
      values = map(GK.vertices_structure(G)) do vertex
        j1, j2 = vertex.first.flag[1], vertex.second.flag[1]
        prod((variables[i] - variables[j1] for i in 1:d1); init=one(first(variables))) *
        prod(
          (variables[3 + i] - variables[3 + j2] for i in 1:d2); init=one(first(variables))
        )
      end
      GK.polynomial_class(G, values)
    else
      throw(ArgumentError("the basis lift must be billey or schur"))
    end
  end
  primitive_edges = map(((1, 0), (0, 1))) do wanted
    first(
      e for e in GK.edges(G)
      if (as_rational(GK.integrate(H1, e)), as_rational(GK.integrate(H2, e))) == wanted
    )
  end
  basis, [GK.curve_class(G, edge) for edge in primitive_edges], collect(primitive_edges)
end

"""Prepare independent native models in a common ambient quotient basis and curve coordinates."""
function prepare_fourfold(name::String; maximum_degree::Int=3, basis_lift::Symbol=:billey)
  Z, V = fourfold_models(name)
  PF.dimension(Z) == 4 || error("the target is not a fourfold")
  X, G = PF.ambient_variety(Z), GK.baseof(V)
  PF.dimension(X) - OS.rank(V) == 4 || error("the native models have different dimensions")
  classical = GW.small_quantum_multiplication(Z; betas=[zeros(Int, PF.picard_rank(X))])
  basis, curves, primitive_edges = if name == "K3-29"
    product_basis(G, V, classical.basis; basis_lift)
  else
    homogeneous_basis(X, G, classical.basis; basis_lift)
  end
  # c5/b8 have degree-five terms; the other targets use their full automatic
  # degree list. A cap is explicit metadata, never presented as a full result.
  index, betas = GW.fano_curve_classes(Z)
  # The Fano index is the gcd, not the vector of degrees on primitive curves.
  # Compute the latter independently in the native GKMtools graph, also for products.
  divisor = GK.first_chern_class(G) - GK.first_chern_class(V)
  anticanonical_degrees = [
    Int(as_rational(GK.integrate(divisor, edge))) for edge in primitive_edges
  ]
  all(>(0), anticanonical_degrees) ||
    error("the native first Chern degrees are not positive")
  foldl(gcd, anticanonical_degrees; init=0) == index || error("the Fano indices disagree")
  name in ("c5", "b8") && filter!(beta -> sum(beta) <= maximum_degree, betas)
  gkm_betas = [
    sum((beta[i] * curves[i] for i in eachindex(curves)); init=zero(parent(first(curves))))
    for beta in betas
  ]
  euler = GK.chern_class(V, OS.rank(V))
  gkm_metric = GK.twisted_intersection_matrix(basis, euler; equivariant=false)
  matrix = [
    as_rational(gkm_metric[i, j]) for
    i in axes(classical.metric, 1), j in axes(classical.metric, 2)
  ]
  matrix == classical.metric || error("the transported classical metrics disagree")
  (; name, Z, V, basis, betas, gkm_betas, gkm_metric,
    codimensions=length.(classical.basis), anticanonical_degrees)
end

"""Compute all requested coefficient matrices using the public GKMtools unmarked API."""
function gkm_operator(context)
  G = GK.baseof(context.V)
  divisor = GK.first_chern_class(G) - GK.first_chern_class(context.V)
  [
    GK.twisted_matrix_small_quantum_product(context.V, divisor;
      basis=context.basis, beta, equivariant=false, show_bar=false)
    for beta in context.gkm_betas
  ]
end

"""Compute the same finite Novikov truncation through the public package API."""
flag_operator(context) = GW.small_quantum_multiplication(context.Z; betas=context.betas)

"""Extract coefficient matrices without folding distinct curve classes together at q=1."""
function flag_coefficients(context, result)
  [map(p -> get(p.terms, Tuple(beta), zero(QQ)), result.matrix) for beta in context.betas]
end

"""Convert output only after timing, preserving the native Oscar arithmetic in the measurement."""
gkm_coefficients(matrices) =
  [[as_rational(M[i, j]) for i in 1:size(M, 1), j in 1:size(M, 2)] for M in matrices]
