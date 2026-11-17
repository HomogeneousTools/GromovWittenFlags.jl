# Conventions and scope

## Schubert classes and curves

The Schubert basis consists of minimal representatives of the right Weyl-group
cosets used by PartialFlagVarieties and Semisimple. Coxeter length is complex
codimension. The first class is the identity and [`point_class`](@ref) has top
codimension. [`schubert_class`](@ref) projects a word to the relevant coset.

Curve coordinates are nonnegative integers in `marked_nodes(X)` order, or a
full ambient Dynkin-node vector that vanishes at unmarked nodes. Thus the
line class on `Gr(2, 4)` may be written `[1]` or `[0, 1, 0]`.
The same ambient coordinates label invariants of a zero locus; a single
ambient class can combine several intrinsic curve classes.

## Dimension and constant maps

For positive curve degree, the required sum of complex codimensions is
`dimension(Z) + c1(TZ)⋅beta + n - 3`. The implementation returns zero when
the dimension axiom fails.

At curve degree zero, the three-point invariant is the classical intersection
product. For fewer than three insertions the API extends this by the same
intersection integral, in particular providing the two-point pairing.
Degree-zero primary invariants with more than three insertions vanish: the
insertions do not provide a top-degree class on the moduli space of marked
curves. No descendant insertions are implemented.

## Defining bundles

Use a `CompletelyReducibleBundle` on the same ambient variety. The engine
expands Levi representations into their weights and requires each weight to
pair nonnegatively with every nonparabolic positive coroot. It then uses the
Euler twist. The weights follow
PartialFlagVarieties' Borel–Weil–Bott convention; the fiber sign is applied
during localization. The defining section must be regular with smooth zero
locus. The package does not test these geometric assumptions.

This weight condition is sufficient for the implemented formulas. Failure of
the condition does not prove that the bundle is nonconvex; the test does not
classify all bundles to which quantum Lefschetz could apply. A negative
weight–coroot pairing is rejected explicitly. A zero torus weight is also
rejected: its Euler factors require a fiber-equivariant regularization that
this package does not implement. In particular, a trivial summand is unsupported.

An actual convex bundle with nonsplit graded pieces can instead be supplied as
[`ConvexEulerTwist`](@ref). The full signed index is assembled and identical
characters are cancelled before numerical specialization. A geometric convexity
certificate and a regular section are caller-supplied hypotheses. See
[Certified computation](@ref) for the contract and an executable model.

## Ambient quantum multiplication

For a zero locus, the ambient intersection pairing may be degenerate.
The implementation chooses independent Schubert restrictions and computes
the projected operator on their quotient by the radical. The `basis` and
`metric` fields describe this choice. Primitive classes and additional
nonambient classes are not added. The ambient subspace need not be preserved
by the full quantum product; this projected matrix alone does not determine
the action on the missing cohomology.

Automatic curve enumeration uses the ambient first Chern coordinates.
If all are positive, it includes every nonnegative vector with first Chern
degree at most `dimension(Z) + 1`, including zero. The returned `index` is
their gcd; it need not equal the intrinsic Fano index after restriction.

If all coordinates vanish, multiplication by the first Chern class is zero.
If some coordinates vanish and others are positive, the numerical bound
cannot give a finite list: supply explicit curve classes. Negative Chern
coordinates also require explicit classes. A supplied list always describes
a truncation, whose completeness is the caller's responsibility.

## Exact arithmetic and resources

Answers lie in `Rational{BigInt}`. The engine retries accidental poles and
requires agreement at two deterministic nonsingular torus points.
Agreement is an implementation check, not a proof of the localization formula.
No floating point or Python call participates. The default backend sums exact
rationals; the opt-in certified backend uses proved bounds and modular
reconstruction, as described in [Certified computation](@ref).
Only a localization pole triggers a retry; disagreement and unrelated errors
propagate. Machine-integer overflow in curve and dimension bounds raises an
error rather than changing the dimension test silently.

Returned Schubert words and Novikov dictionaries have their own storage.
Treat a polynomial's `terms` dictionary as read-only; use arithmetic to form
new polynomials. Private graph callbacks instead borrow reusable buffers,
as described in [Internals](@ref).

The number of decorated trees grows quickly with the sum of curve coordinates,
and Weyl coset enumeration can be large in high rank. Structural caches can
be released with [`clear_localization_caches!`](@ref). Evaluation-specific
caches are local and become collectible after a computation.
