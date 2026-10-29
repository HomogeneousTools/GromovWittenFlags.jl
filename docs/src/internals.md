# Internals

This page describes the implementation contracts used by contributors.
Underscored bindings are private and may change between releases.
The mathematical reference is the root-theoretic localization construction in
[arXiv:2609.21914v2](https://arxiv.org/abs/2609.21914v2).
The geometric assumptions and the ambient projection are described in
[Conventions and scope](@ref).

## Following a computation

The public boundary validates geometry, converts curve coordinates, and resolves
insertions to integer Schubert indices. Arithmetic then runs on concrete cached
root data. The files follow this order:

| Layer | Files | Responsibility |
|:--|:--|:--|
| Geometry | `Schubert.jl`, `Bundles.jl` | Cosets, acted roots, curves, and defining weights |
| Fixed graphs | `FreeTrees.jl`, `TreeSymmetries.jl`, `Colorings.jl`, `BudgetColorings.jl`, `DecoratedTrees.jl` | Tree orbits, stabilizers, and borrowed workspaces |
| Numeric kernels | `ExactSums.jl`, `Localization.jl` | Exact torus values, restrictions, and vertex/edge contributions |
| Shared insertions | `QuantumBlocks.jl` | Two-point blocks with shared graph enumeration |
| Public invariants | `Invariants.jl` | Input guards, dimension axiom, and specialization orchestration |
| Quantum assembly | `Novikov.jl`, `ExactLinearAlgebra.jl`, `Quantum.jl` | Ambient quotient, metric solves, and sparse coefficients |

The low-level localization file does not call the higher-level pair-block
dispatcher. [`GromovWittenFlags._checked_localization`](@ref) chooses that
optimization in the invariant layer. Restriction tables belong to the numeric
layer because both general invariants and pair blocks use them.

Several integers have different meanings. Preserve these distinctions when
adding a kernel; names such as `root`, `class`, and `index` alone are insufficient
at a boundary between these conventions:

| Quantity | Meaning |
|:--|:--|
| Fixed-point color / full Schubert index | Position in `data.representatives` |
| Local root index `q` | Position in `data.nonparabolic_indices` and `data.curve_classes` |
| Ambient positive-root index | `data.nonparabolic_indices[q]` |
| Block class index `i` | Position in a supplied `classes` map; `classes[i]` is a full Schubert index |
| Curve coordinate `beta[p]` | Degree at ambient Dynkin node `marked_nodes(X)[p]` |
| Cover degree | Positive multiplicity on a single invariant edge |

## Structural data and ownership

[`GromovWittenFlags._GWData`](@ref) belongs to one marked Dynkin type.
It caches minimal right-coset representatives, inverse representatives,
codimensions, acted roots, reduced-word roots, curve classes, root–coroot
pairings, and the moment graph. Fixed-point colors are indices into its
representative vector. An edge root is an index into `nonparabolic_indices`,
expressed in the source color's root frame.

The heterogeneous global variety cache ends at a typed function barrier.
Arithmetic kernels dispatch on the concrete root datum and scalar type.
Torus-independent root coordinates use fixed-size integer vectors;
specialization evaluates their integral linear forms before converting to
exact rationals. Positive-degree quantum pair blocks use Nemo's FLINT-backed
`QQFieldElem` internally, then convert checked coefficients to `Rational{BigInt}`
before matrix assembly. Classical and scalar-invariant localization continue
to use `Rational{BigInt}`. Machine-integer degree bounds use checked arithmetic.

The heterogeneous variety dictionary and large-tree tuple encodings are
deliberate boundaries. Hot arithmetic does not use `Any`-typed fields. A torus
evaluation stores the fiber values and splittings needed by its own target;
it does not retain the input weight list or an unused copy of tangent Euler
classes. Twisted and untwisted evaluations share the same concrete scalar-based
layout.

Public Schubert constructors copy the Weyl word. Mutating a returned word
therefore cannot corrupt the cache. A public [`NovikovPolynomial`](@ref)
constructor likewise owns a normalized copy of its input dictionary.
Internal polynomial assembly mutates only freshly allocated result entries.

Tree shapes and automorphism groups are also structural caches.
Their arrays are shared and treated as read-only. Cache insertion is protected
by locks; each computation owns its arithmetic caches and mutable workspace.
[`clear_localization_caches!`](@ref) releases the structural cache entries.
Existing computations retain their own references to data already obtained.

Mutating workspace helpers use `!` and document which buffers they overwrite.
Memoized readers such as `_billey` and `_edge_factor` may populate their private
caches; their mathematical values and geometry remain unchanged. These effects
are explicit in their docstrings. Callback consumers must treat borrowed graph
arrays as read-only, as well as avoiding references that outlive a callback.

## Schubert restrictions

[`GromovWittenFlags._billey`](@ref) evaluates the reduced-subword formula.
For a target class, [`GromovWittenFlags._billey_plan`](@ref) builds an automaton
of all possible reduced prefixes of that target. It obtains these by repeated
right descents; the dependency calls this operation `bruhat_descendants`,
although these immediate descendants are the right weak-order descents.

Each state records the next state for a simple generator, or zero if that
transition cannot reach the target. At each letter of the fixed-point word,
the recurrence updates states in descending length. This order prevents using
one letter twice. The identity restriction is one; shorter fixed-point words
give zero, and the diagonal restriction is the product of the word roots.
Restrictions are memoized at each torus point and shared by all insertion
pairs in a block.

The tests compare this recurrence against literal subsets of reduced-word
positions, including non-simply-laced root systems and products.

## Tree enumeration and symmetry

Free trees use the inherited Wright–Richmond–Odlyzko–McKay successor algorithm.
Each yielded depth sequence has its own storage. A cached tree shape records
edges, incidence, and a parent-before-child traversal.

Canonical colorings enumerate graph homomorphisms from the tree to the
moment graph modulo uncolored tree automorphisms. The iterator keeps its
successor state separate from the owned vectors it yields. When a color
advances, complete subtrees in the changed suffix reset in place, using cached
parent indices; no subtree level slices or temporary color vectors are needed.
Rooted enumeration fixes vertex one's color and allows only root-preserving
symmetries. Tests compare both rooted and unrooted orbits with brute force.

Decorated-tree enumeration uses a preorder callback search with a remaining
curve budget. Each assigned edge spends its all-unit curve class; if any
coordinate exceeds the budget, every positive-cover completion is impossible.
Identical-sibling prefix comparisons and the unmarked central flip preserve
the same color orbits as the original iterator, which remains a test oracle.
The callback borrows the color/root buffers and minimum class. Cover degrees
are then enumerated modulo the subgroup preserving the colors. For small shapes,
permutation groups give degree-orbit representatives and stabilizer orders.
For larger shapes, recursive canonical encodings replace explicit groups.
The stabilizer order contributes the decorated-tree automorphism denominator;
the separate cyclic cover denominator is included in each edge factor.

The cover-degree search streams positive integer assignments.
Suffix sums reserve the minimum class needed by the remaining edges, and the
last edge degree is solved directly. The all-unit case bypasses degree-orbit
keys. A rooted calculation restricts automorphisms to those fixing the
distinguished vertex and uses one rooted shape per uncolored vertex orbit.

The cover and finite-curve searches use named recursive workers with explicit
scratch arguments. A recursive local function can capture its own binding in a
`Core.Box`; passing state explicitly avoids that dispatch and allocation cost.
The degree worker restores its remaining budget after each branch and copies
nothing for a callback. The finite-curve worker copies a vector only when it
appends a complete result. The large-tree encoder returns its key and symmetry
order together, so it needs no shared mutable counter.

[`GromovWittenFlags._foreach_unmarked_decorated_tree`](@ref) and
[`GromovWittenFlags._foreach_rooted_decorated_tree`](@ref) are synchronous
callback APIs. **All workspace arrays passed to a callback are borrowed.**
Copy any array that must survive the callback. Structural arrays must never
be mutated. Exceptions propagate on the calling task, including when a pole
aborts a specialization; no suspended producer task retains the workspace.

## Localization factors and labeled marks

[`GromovWittenFlags._GWEvaluation`](@ref) binds one exact torus point to one
defining bundle. Fundamental-weight coordinates are converted to simple-root
coordinates before evaluation. The minus sign in the fiber character converts
PartialFlagVarieties' Borel–Weil–Bott convention to the actual fiber action.
Root–coroot pairings and weight pairings determine cover splittings.

At a fixed point, `vertex_euler` is the tangent Euler class divided by the
bundle Euler class. The endpoint factors from all edges and the vertex
denominator combine into `vertex_euler^(valence - 1)`. Collecting them first
avoids multiplying large rational factors that immediately cancel.
[`GromovWittenFlags._edge_factor`](@ref) contains the moving-edge and bundle
factors and the cyclic cover denominator, but no endpoint Euler classes.

For incident tangent weights `w_i`, put `s = sum(1 / w_i)`.
The zero-mark descendant factor is `prod(1 / w_i) * s^(valence - 3)`;
each labeled mark at the vertex multiplies it by `s` and the insertion's
Schubert restriction. Leaf and bivalent factors are cached by fixed-point
color and unordered local root/cover pairs. Higher valences are computed
directly to keep the cache size under control.

For one insertion, summing its possible vertices gives
`sum(mark_weight[v] * restriction[color[v]])`.
Labeled marks choose their vertices independently, so their sums multiply.
Orbit–stabilizer identifies the sum over all labeled placements divided by
the unmarked automorphism order with the sum over marked orbits divided by
their marked stabilizer orders. This removes Cartesian mark enumeration;
repeated insertions reuse one vertex sum. A two-point block shares these sums
across all requested pairs and across the two torus evaluations.
When a singleton pair has small fixed-point support, the rooted algorithm
anchors one mark and sums the second mark over the vertices.

Degree-zero integration has a separate fixed-point sum and the classical
pairing convention. It never passes a valence-zero graph to the positive-degree
vertex formula. Tests compare the mark factorization with explicit labeled
placements and marked stabilizers on small graphs.

## Poles and exact comparisons

Julia rational division can produce an infinite rational such as `1//0`.
[`GromovWittenFlags._divide_or_pole`](@ref) explicitly tests denominators;
the final scalar or matrix is also checked for finiteness.
All feasible covered-edge factors are preflighted before restriction-support
pruning. Thus skipping a zero insertion contribution cannot hide an edge pole.
Vertex factors are evaluated before pruning their mark sums.

Two deterministic nonsingular torus specializations must agree exactly.
The denominator and finiteness guards signal an explicit pole with the private
[`GromovWittenFlags._LocalizationPole`](@ref) exception. Only this signal retries
the whole block at the next pair. Eight small-coordinate pairs are tried first,
followed by all eight original larger-coordinate pairs. Ordinary `DivideError`s
propagate immediately, just as domain errors and other programming exceptions do.
Disagreement is an error and never retries. Agreement checks consistency
of the implementation; it is not a proof of the localization formula.

Graph contributions use [`GromovWittenFlags._ExactSum`](@ref), a streaming
balanced sum. Each occupied binary slot holds a block of equally many nonzero
terms; adding a term merges equal-size blocks like a binary carry. This keeps
most rational gcd operations between small partial sums instead of repeatedly
combining one graph with the large denominator of the entire prefix. The
accumulator retains logarithmically many blocks and changes no coefficients.

## Quantum assembly and exact algebra

The classical metric and first-Chern pairing share one restriction table and
one target-weight expansion. Codimension gates skip entries known to vanish.
Incremental rational row elimination selects the first independent ambient
classes. For a symmetric metric these classes give a nonsingular principal
submatrix, representing the quotient by its radical.

One exact LU factorization of this reduced metric is reused for the classical
and every positive-degree coefficient matrix. The divisor axiom converts a
first-Chern insertion to its curve degree times a two-point block. Solving
against the metric raises the index; sparse Novikov entries receive only
nonzero coefficients. This avoids forming an explicit inverse or repeatedly
copying entire polynomial matrices.

Automatic finite curve enumeration requires positive first-Chern coordinates.
Mixed zero directions and negative directions need explicit truncations.
No primitive cohomology is added to this ambient calculation.

## Private reference

The following docstrings are checked by Documenter. Their contracts explain
arguments, mathematical conventions, and storage lifetimes where applicable.

```@autodocs
Modules = [GromovWittenFlags]
Public = false
Private = true
Order = [:type, :constant, :function]
```
