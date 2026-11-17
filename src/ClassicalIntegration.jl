"""
    equivariant_flag_integral(integrand, X; bundles=(), classes=(), torus_seeds=(0,3))

Integrate a top-degree equivariant characteristic expression on a flag variety,
using exact rational localization at two torus points. `integrand(fibers, restrictions)`
receives the lists of actual fiber characters of each supplied bundle and the
restrictions of the supplied Schubert classes at one fixed point. It may return
a rational scalar or an array of rational scalars. Its output must be homogeneous
of complex degree `dimension(X)`, including any Euler class of a zero locus.
The callback is responsible for this degree condition. The two answers must agree.
Zero fiber weights are supported; no Euler division by the supplied bundles occurs.
"""
function equivariant_flag_integral(integrand, X::PartialFlagVariety;
  bundles=(), classes=(), torus_seeds=(0, 3))
  seeds = Int.(collect(torus_seeds))
  length(seeds) == 2 && length(unique(seeds)) == 2 && all(s -> 0 <= s < 32, seeds) ||
    throw(ArgumentError("choose two distinct torus seeds in 0:31"))
  data = _gw_data(X)
  all(E -> E isa CompletelyReducibleBundle && variety(E) == X, bundles) ||
    throw(ArgumentError("all bundles must belong to the supplied flag variety"))
  weights = [_bundle_weights(E) for E in bundles]
  indices = [_class_index(data, s) for s in classes]
  Cinv = _BigRational.(cartan_matrix_inverse(dynkin_type(X)))
  answers = map(seeds) do seed
    evaluation = _gw_evaluation(data, seed)
    table = _restriction_table(evaluation, indices)
    total = nothing
    for color in eachindex(data.representatives)
      w_inverse = data.inverse_representatives[color]
      fibers = [
        [-sum((Cinv * coefficients(w * w_inverse)) .* evaluation.point) for w in ws]
        for ws in weights
      ]
      value = integrand(fibers, @view(table[color, :])) / evaluation.vertex_euler[color]
      total = total === nothing ? value : total + value
    end
    total
  end
  answers[1] == answers[2] || error("classical integrals disagree at the two torus points")
  first(answers)
end
