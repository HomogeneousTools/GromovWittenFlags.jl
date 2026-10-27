"""
Prepare the classical data once for a second, deliberately favorable GKMtools comparison.

The public matrix API recomputes its symbolic pairing for each curve class.
This adapter retains its own independently computed classical matrix and metric.
All positive-degree integration still uses the public GKMtools evaluator.
Its preparation cost is recorded separately, rather than hidden in a speedup.
"""
function prepare_gkm_operator(context)
  G = GK.baseof(context.V)
  divisor = GK.first_chern_class(G) - GK.first_chern_class(context.V)
  classical = GK.twisted_matrix_small_quantum_product(context.V, divisor;
    basis=context.basis, beta=zero(parent(first(context.gkm_betas))),
    equivariant=false, show_bar=false)
  (; divisor, classical, inverse_metric=OS.inv(context.gkm_metric),
    twist=GK.virtual_zero_section(context.V))
end

"""
Assemble coefficients with a retained classical pairing and public unmarked integration.

The dimension constraint for a fourfold is 1+codim(a)+codim(b)=4+c₁⋅β.
Only those pairs are submitted; the upstream matrix API applies the same test.
Every retained coefficient is checked against the independent package result.
This is a benchmark adapter, not a proposed production backend or a change to GKMtools.
"""
function gkm_prepared_operator(context)
  prepared = context.prepared
  G = GK.baseof(context.V)
  n = length(context.basis)
  map(zip(context.betas, context.gkm_betas)) do (coordinates, beta)
    iszero(beta) && return prepared.classical
    degree = sum(context.anticanonical_degrees .* coordinates)
    pairs = [
      (j, k) for j in 1:n for k in j:n
      if 1 + context.codimensions[j] + context.codimensions[k] == 4 + degree
    ]
    matrix = OS.zero_matrix(OS.QQ, n, n)
    isempty(pairs) && return matrix
    products = [[prepared.divisor, context.basis[j], context.basis[k]] for (j, k) in pairs]
    invariants = GK.gromov_witten_nomarks(G, beta, products, 0, prepared.twist;
      fast_mode=true, show_bar=false, known_zero=fill(false, length(pairs)))
    for ((j, k), value) in zip(pairs, invariants)
      matrix[j, k] = matrix[k, j] = value
    end
    prepared.inverse_metric * matrix
  end
end
