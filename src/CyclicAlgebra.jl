"""Evaluate a rational polynomial at an exact matrix by Horner's rule."""
function _cert_matrix_polynomial(p, A)
  n = Nemo.nrows(A)
  identity = Nemo.identity_matrix(Nemo.QQ, n)
  result = Nemo.zero_matrix(Nemo.QQ, n, n)
  for k in reverse(0:Nemo.degree(p))
    result = result * A + Nemo.coeff(p, k) * identity
  end
  result
end

"""
    cyclic_quantum_algebra(operator, unit; metric=nothing)

Recover an exact commutative algebra from multiplication by a cyclic element.
Inputs are a rational matrix at a specified Novikov specialization and its unit
vector. Return the characteristic/minimal polynomial, Krylov determinant, all
multiplication operators, and rational primary projectors. Reject a noncyclic
unit. Check associativity, commutativity, projector identities, and, when supplied,
Frobenius self-adjointness. This computes spectral algebra data, not geometric
Hodge-atom hypotheses or genericity of big quantum deformations.
"""
function cyclic_quantum_algebra(
  operator::AbstractMatrix, unit::AbstractVector; metric=nothing
)
  n, m = size(operator)
  n == m == length(unit) || throw(DimensionMismatch("operator and unit"))
  A = Nemo.matrix(Nemo.QQ, operator)
  e = Nemo.matrix(Nemo.QQ, reshape(unit, n, 1))
  identity = Nemo.identity_matrix(Nemo.QQ, n)
  zero = Nemo.zero_matrix(Nemo.QQ, n, n)
  columns = [A^k * e for k in 0:(n - 1)]
  K = reduce(hcat, columns)
  determinant = Nemo.det(K)
  iszero(determinant) && throw(ArgumentError("the supplied unit is not cyclic"))
  inverse = inv(K)
  R, t = Nemo.polynomial_ring(Nemo.QQ, "t")
  polynomial = Nemo.charpoly(R, A)
  operators = map(1:n) do i
    column = inverse[:, i:i]
    B = zero
    for k in 0:(n - 1)
      B += column[k + 1, 1] * A^k
    end
    B * e == identity[:, i:i] || error("cyclic unit check failed")
    B
  end
  for i in 1:n, j in 1:n
    operators[i][:, j:j] == operators[j][:, i:i] || error("commutativity check failed")
    operators[i] * operators[j] == sum(
      (operators[i][k, j] * operators[k] for k in 1:n); init=zero) ||
      error("associativity check failed")
  end
  if metric !== nothing
    G = Nemo.matrix(Nemo.QQ, metric)
    all(B -> transpose(B)*G == G*B, operators) || error("Frobenius pairing check failed")
  end
  blocks = Any[]
  for (factor, exponent) in Nemo.factor(polynomial)
    primary = factor^exponent
    quotient = Nemo.divexact(polynomial, primary)
    inverse_polynomial = Nemo.invmod(quotient, primary)
    projector_polynomial = rem(quotient * inverse_polynomial, polynomial)
    projector = _cert_matrix_polynomial(projector_polynomial, A)
    projector^2 == projector || error("invalid primary projector")
    Nemo.rank(projector) == Nemo.degree(primary) || error("invalid primary rank")
    push!(
      blocks,
      (; factor, exponent, dimension=Nemo.degree(primary),
        projector_polynomial, projector, local_unit=projector*e),
    )
  end
  sum((b.projector for b in blocks); init=zero) == identity ||
    error("incomplete projectors")
  for i in eachindex(blocks), j in 1:(i - 1)
    blocks[i].projector * blocks[j].projector == zero ||
      error("nonorthogonal primary projectors")
  end
  (; characteristic_polynomial=polynomial, minimal_polynomial=polynomial,
    krylov_determinant=determinant, multiplication_operators=operators, blocks)
end
