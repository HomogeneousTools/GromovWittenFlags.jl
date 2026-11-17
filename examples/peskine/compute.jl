module PeskineExample
# No saved basis, residues, bounds, quantum coefficients, or Python inputs are read.
const PESKINE_PROCESS_STARTED = time_ns()
using PartialFlagVarieties, GromovWittenFlags
import Nemo, Semisimple, LinearAlgebra, TOML, SHA

const PESKINE_CODIMS = [0, 1, 2, 2, 3, 4]

function peskine_model()
  F = flag_variety(10, [1, 4])
  L, B = tautological_bundles(F)
  C = CompletelyReducibleBundle(F, Semisimple.fundamental_weight(dynkin_type(F), 9))
  graded =
    tensor_product(dual(L), exterior_power(dual(B), 2) +
                            tensor_product(dual(B), dual(C))) + 2dual(L)
  twist = ConvexEulerTwist(
    graded;
    certificate=
    "The actual rank-21 bundle is Ldual tensor (wedge2(V/L)dual / wedge2 Cdual). " *
    "Contraction from wedge3 Vdual makes it globally generated. " *
    "Two further Ldual sections give a smooth regular Peskine fourfold.",
  )
  F, L, B, twist
end

function peskine_preferred_basis(F, L, B, twist, result)
  codims = [length(Semisimple.word(s)) for s in result.basis]
  cross = equivariant_flag_integral(F;
    bundles=(twist.graded, L, B), classes=result.basis) do fibers, restrictions
    h = -sum(fibers[2])
    b = (sum(fibers[3])^2 - sum(w^2 for w in fibers[3])) / 2
    tautological = [one(h), h, h^2, b, h^3, h^4]
    euler = prod(fibers[1])
    [
      if codims[i] + PESKINE_CODIMS[j] == 4
        euler * restrictions[i] * tautological[j]
      else
        zero(h)
      end for i in 1:6, j in 1:6
    ]
  end
  change = result.metric \ cross
  metric = transpose(change) * result.metric * change
  iszero(LinearAlgebra.det(change)) && error("preferred basis is dependent")
  change, metric
end

function peskine_full_trace()
  F = flag_variety(10, [1, 2, 7])
  L, M, B = tautological_bundles(F)
  U = L + M
  C = CompletelyReducibleBundle(F, Semisimple.fundamental_weight(dynkin_type(F), 9))
  H = dual(L)
  defining =
    tensor_product(exterior_power(dual(U), 2), dual(B)) +
    tensor_product(dual(U), exterior_power(dual(B), 2)) + 2dual(U)
  H3 = tensor_product(tensor_product(H, H), H)
  positive = 8H + tensor_product(dual(C), tensor_product(H3, H)) + tensor_product(L, H3)
  negative = tensor_product(L, H) + tensor_product(U+B, H3) + tensor_product(H, M)
  @assert dimension(F) == 32 && PartialFlagVarieties.rank(defining) == 29
  equivariant_flag_integral(F; bundles=(defining, positive, negative)) do fibers, _
    powers = [sum(w^j for w in fibers[2]) - sum(w^j for w in fibers[3]) for j in 1:3]
    prod(fibers[1]) * (powers[1]^3 - 3powers[1]*powers[2] + 2powers[3]) / 6
  end
end

function peskine_line_moments()
  F = flag_variety(10, [2, 7])
  U, B = tautological_bundles(F)
  C = CompletelyReducibleBundle(F, Semisimple.fundamental_weight(dynkin_type(F), 9))
  defining =
    tensor_product(exterior_power(dual(U), 2), dual(B)) +
    tensor_product(dual(U), exterior_power(dual(B), 2)) + 2dual(U)
  equivariant_flag_integral(F; bundles=(defining, U, C)) do fibers, _
    a, c = -sum(fibers[2]), sum(fibers[3])
    u2 = (a^2-sum(w^2 for w in fibers[2]))/2
    c2 = (c^2-sum(w^2 for w in fibers[3]))/2
    prod(fibers[1]) .* [a^2, u2, a*c, c^2, c2]
  end
end

peskine_strings(A) = [[string(A[i, j]) for j in axes(A, 2)] for i in axes(A, 1)]

function peskine_main(directory, workers; maximum_degree=5)
  mkpath(directory)
  F, L, B, twist = peskine_model()
  progress = event -> begin
    println(event)
    flush(stdout)
  end
  started = time_ns()
  result = small_quantum_multiplication(F; bundle=twist,
    betas=[[d, 3d] for d in 0:maximum_degree], arithmetic=:certified,
    workers, worker_heap_size="4G", checkpoint=joinpath(directory, "checkpoints"), progress,
  )
  quantum_seconds = (time_ns()-started)/1e9
  change, metric = peskine_preferred_basis(F, L, B, twist, result)
  inverse = inv(change)
  blocks = Matrix{Rational{BigInt}}[]
  invariants = Any[]
  for d in 0:maximum_degree
    raw = [get(p.terms, (d, 3d), big(0)//1) for p in result.matrix]
    block = inverse * raw * change
    transpose(block)*metric == metric*block || error("Frobenius pairing mismatch")
    all(
      iszero(block[i, j]) || PESKINE_CODIMS[j]+1-PESKINE_CODIMS[i] == d
      for i in 1:6, j in 1:6
    ) || error("quantum grading mismatch")
    push!(blocks, block)
    if d > 0
      pairings = metric*block/d
      for i in 2:6, j in i:6
        PESKINE_CODIMS[i]+PESKINE_CODIMS[j] == 3+d || continue
        push!(
          invariants, Dict("degree"=>d, "i"=>i, "j"=>j, "value"=>string(pairings[i, j]))
        )
      end
    end
  end
  report = Dict{String,Any}("status"=>maximum_degree == 5 ? "complete" : "partial",
    "julia_version"=>string(VERSION), "workers"=>workers,
    "basis"=>["1", "h", "h^2", "c2(B3)", "h^3", "h^4"],
    "codimensions"=>PESKINE_CODIMS, "metric"=>peskine_strings(metric),
    "classical_basis_change"=>peskine_strings(change),
    "schubert_basis_words"=>[collect(Semisimple.word(s)) for s in result.basis],
    "coefficient_blocks"=>[
      Dict("degree"=>d, "matrix"=>peskine_strings(blocks[d + 1])) for d in 0:maximum_degree
    ],
    "two_point_invariants"=>invariants, "certification"=>result.certification,
    "quantum_seconds"=>quantum_seconds,
    "inputs"=>"only the flag/bundle model, curve degrees, and stated geometric Hodge ranks; no previous computation files",
  )
  if maximum_degree == 5
    operator = sum(blocks)
    algebra = cyclic_quantum_algebra(operator, [1, 0, 0, 0, 0, 0]; metric)
    trace = peskine_full_trace()
    moments = peskine_line_moments()
    trace == 5moments[1]-7moments[3]+8moments[4]-16moments[5] ||
      error("two trace calculations disagree")
    roots = [
      -Nemo.coeff(block.factor, 0)/Nemo.coeff(block.factor, 1) for block in algebra.blocks
    ]
    all(Nemo.degree(block.factor)==1 for block in algebra.blocks) ||
      error("nonrational spectrum")
    solutions = [
      (a, b) for a in roots, b in roots if
                 LinearAlgebra.tr(operator)+22a+32b == trace
    ]
    length(solutions)==1 ||
      error("trace does not assign the two primitive actions uniquely")
    alpha, beta = only(solutions)
    R = Nemo.parent(algebra.characteristic_polynomial)
    t = Nemo.gen(R)
    full = algebra.characteristic_polynomial*(t-alpha)^22*(t-beta)^32
    report["characteristic_polynomial_at_q_1"] = string(algebra.characteristic_polynomial)
    report["minimal_polynomial_at_q_1"] = string(algebra.minimal_polynomial)
    report["krylov_determinant_at_q_1"] = string(algebra.krylov_determinant)
    report["multiplication_operators_at_q_1"] = peskine_strings.(
      algebra.multiplication_operators
    )
    report["full_degree_one_trace"] = string(trace)
    report["line_surface_moments"] = string.(moments)
    report["inherited_scalar"] = string(alpha)
    report["vanishing_scalar"] = string(beta)
    report["full_characteristic_polynomial_at_q_1"] = string(full)
    report["spectral_blocks"] = [
      Dict(
        "eigenvalue_at_q_1"=>string(roots[i]), "ambient_dimension"=>block.dimension,
        "projector_at_q_1"=>peskine_strings(block.projector),
        "local_unit_at_q_1"=>peskine_strings(block.local_unit),
        "hodge_classes"=>block.dimension,
        "inherited_rank"=>roots[i]==alpha ? 22 : 0,
        "vanishing_rank"=>roots[i]==beta ? 32 : 0,
        "full_dimension"=>block.dimension+(roots[i]==alpha ? 22 : 0)+(
          roots[i]==beta ? 32 : 0
        ),
      ) for (i, block) in enumerate(algebra.blocks)
    ]
    report["atom_scope"] = "Very-general coarse Hodge atoms, using absolute irreducibility, Hodge-equivariance, and the regular F-manifold genericity argument. These geometric hypotheses are not inferred numerically."
  end
  report["process_elapsed_seconds"] = (time_ns()-PESKINE_PROCESS_STARTED)/1e9
  report["driver_sha256"] = bytes2hex(SHA.sha256(read(@__FILE__)))
  open(joinpath(directory, "peskine.toml"), "w") do io
    TOML.print(io, report; sorted=true)
  end
  println("RESULT ", joinpath(directory, "peskine.toml"))
  println("QUANTUM_SECONDS ", quantum_seconds)
  println("PROCESS_ELAPSED_SECONDS ", report["process_elapsed_seconds"])
  report
end

if abspath(PROGRAM_FILE)==@__FILE__
  length(ARGS)>=1 || error("usage: compute.jl OUTPUT_DIRECTORY [WORKERS=4] [MAX_DEGREE=5]")
  peskine_main(abspath(ARGS[1]), length(ARGS)>1 ? parse(Int, ARGS[2]) : 4;
    maximum_degree=length(ARGS)>2 ? parse(Int, ARGS[3]) : 5)
end

end
