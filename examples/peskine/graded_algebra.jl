# Finite exact export from a completed computation; no localization is repeated.
using GromovWittenFlags
import TOML, SHA

function graded_rational(s)
  parts = split(replace(s, "//"=>"/"), "/")
  parse(BigInt, parts[1]) // (length(parts)==1 ? big(1) : parse(BigInt, parts[2]))
end

function graded_quantum_operators(specialized, codimensions, metric)
  P = NovikovPolynomial{1}
  n = length(codimensions)
  operators = map(1:n) do i
    coefficients = specialized[i]
    [
      begin
        coefficient = coefficients[k, j]
        exponent = codimensions[i]+codimensions[j]-codimensions[k]
        iszero(coefficient) || exponent >= 0 || error("negative Novikov exponent")
        iszero(coefficient) ? zero(P) : P(Dict((exponent,)=>coefficient))
      end for k in 1:n, j in 1:n
    ]
  end
  identity = [P(i==j ? 1 : 0) for i in 1:n, j in 1:n]
  operators[1]==identity || error("graded unit check failed")
  for i in 1:n
    permutedims(operators[i])*metric == metric*operators[i] ||
      error("graded Frobenius check failed")
    for j in 1:n
      operators[i][:, j] == operators[j][:, i] || error("graded commutativity check failed")
      product = [
        sum((operators[i][k, j]*operators[k][row, column] for k in 1:n);
          init=zero(P)) for row in 1:n, column in 1:n
      ]
      operators[i]*operators[j] == product || error("graded associativity check failed")
    end
  end
  operators
end

function export_graded_algebra(directory)
  path = joinpath(directory, "peskine.toml")
  report = TOML.parsefile(path)
  n = length(report["codimensions"])
  rational_matrix(rows) = [graded_rational(rows[i][j]) for i in 1:n, j in 1:n]
  specialized = rational_matrix.(report["multiplication_operators_at_q_1"])
  metric = rational_matrix(report["metric"])
  operators = graded_quantum_operators(specialized, report["codimensions"], metric)
  for block in report["coefficient_blocks"]
    expected = rational_matrix(block["matrix"])
    d = block["degree"]
    actual = [get(p.terms, (d,), big(0)//1) for p in operators[2]]
    actual==expected || error("exported h multiplication disagrees with certified blocks")
  end
  result = Dict("status"=>"passed", "novikov_variable"=>"q1 of intrinsic degree one",
    "basis"=>report["basis"], "codimensions"=>report["codimensions"],
    "source_peskine_sha256"=>bytes2hex(SHA.sha256(read(path))),
    "checks"=>[
      "unit",
      "nonnegative exponents",
      "Frobenius pairing",
      "commutativity",
      "associativity",
      "certified h multiplication",
    ],
    "multiplication_operators"=>[
      Dict("basis_index"=>i, "matrix"=>[[string(A[k, j]) for j in 1:n] for k in 1:n])
      for (i, A) in enumerate(operators)])
  open(joinpath(directory, "quantum-algebra.toml"), "w") do io
    TOML.print(io, result; sorted=true)
  end
  result
end
