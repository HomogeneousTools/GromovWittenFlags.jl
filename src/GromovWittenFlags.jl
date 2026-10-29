"""
    GromovWittenFlags

Exact genus-zero primary Gromov--Witten invariants of partial flag varieties
and Euler-twisted invariants of supported homogeneous zero loci.
Construct varieties and bundles with PartialFlagVarieties.jl.
"""
module GromovWittenFlags

using Combinatorics: permutations
using LinearAlgebra: Diagonal, lu
using Nemo: QQFieldElem
using PartialFlagVarieties:
  PartialFlagVariety, MarkedDynkinType, CompletelyReducibleBundle, IrrepLevi,
  ZeroLocus, marked_dynkin_type, dynkin_type, marked_nodes, unmarked_nodes,
  dimension, picard_rank, variety, components, det, p_dominant_weight,
  is_borel, semisimple_part, levi_type, levi_permutation,
  ambient_variety, defining_bundle
using Semisimple:
  DynkinType, WeylGroup, WeylGroupElem, WeightLatticeElem, RootSpaceElem,
  RootSystem, rank, weyl_group, right_coset_reps, fundamental_weight,
  coefficients, positive_roots, positive_coroots, reflect,
  cartan_matrix_inverse, freudenthal_formula, simple_root, word, gen,
  bruhat_descendants
using StaticArrays: SVector

export schubert_classes, schubert_class, point_class
export expected_gromov_witten_degree, gromov_witten_invariant
export fano_curve_classes, small_quantum_multiplication
export NovikovPolynomial, evaluate, clear_localization_caches!

"""Exact coefficient field; localization cancellations must never pass through floating point."""
const _BigRational = Rational{BigInt}

include("Schubert.jl")
include("Bundles.jl")
include("FreeTrees.jl")
include("TreeSymmetries.jl")
include("Colorings.jl")
include("BudgetColorings.jl")
include("DecoratedTrees.jl")
include("ExactSums.jl")
include("Localization.jl")
include("QuantumBlocks.jl")
include("Invariants.jl")
include("Novikov.jl")
include("ExactLinearAlgebra.jl")
include("Quantum.jl")
include("Caches.jl")

end
