"""Serial full-matrix comparison with an unchanged archived HPC numerical engine."""
module HPCComparison

using Dates: UTC, now
using LinearAlgebra: BLAS, diag
using SHA: sha256
using TOML
import GromovWittenFlags
import PartialFlagVarieties as PFV
import Semisimple

"""Dependency imports formerly supplied by the parent PartialFlagVarieties module."""
module OldEngine
  using Combinatorics: permutations
  using LinearAlgebra: Diagonal, inv
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
end

"""Flush phase boundaries so long runs remain observable through redirected logs."""
function progress(label, backend, phase)
  println(now(UTC), " ", label, " ", backend, " ", phase)
  flush(stdout)
end

"""Write a checkpoint by renaming a complete TOML file in the same directory."""
function checkpoint(path, value)
  mkpath(dirname(path))
  temporary = path * ".tmp." * string(getpid())
  open(temporary, "w") do io
    TOML.print(io, value; sorted=true)
  end
  mv(temporary, path; force=true)
end

"""Canonical exact scalar representation shared with the archived HPC JSON outputs."""
function rational_text(x::Rational)
  denominator(x) == 1 && return string(numerator(x))
  string(numerator(x), '/', denominator(x))
end

"""Serialize every nonzero Novikov coefficient in lexicographic exponent order."""
function polynomial_text(polynomial)
  terms = sort!(collect(polynomial.terms); by=first)
  join((join(beta, ',') * "=" * rational_text(value) for (beta, value) in terms), ';')
end

"""Label a type-A Schubert representative by its stabilizer-weight coordinates."""
function basis_key(data, w, factors)
  omega = collect(Semisimple.coefficients(data.stabilizer_weight * inv(w)))
  result = Int[]
  offset = 0
  for factor in factors
    n = Int(last(factor))
    block = omega[(offset + 1):(offset + n - 1)]
    epsilon = [sum(block[j:end]; init=0) for j in 1:(n - 1)]
    push!(epsilon, 0)
    append!(result, n .* epsilon .- sum(epsilon))
    offset += n - 1
  end
  join(result, ',')
end

"""Apply one sparse coefficient matrix of the anticanonically specialized connection."""
function apply_coefficient(entries, vector)
  result = zeros(Rational{BigInt}, length(vector))
  for (i, j, value) in entries
    iszero(vector[j]) || (result[i] += value * vector[j])
  end
  result
end

"""
Compute the fourfold period through `order`, using the HPC recurrence exactly.
Grading makes the degree-zero coefficient nilpotent of order at most five;
four additional applications therefore invert the classical part at each degree.
This check is outside the quantum-matrix timing.
"""
function period_coefficients(matrix, codimensions, chern_coordinates, order=100)
  n = size(matrix, 1)
  top = only(findall(==(4), codimensions))
  Q = Rational{BigInt}
  powers = [Tuple{Int,Int,Q}[] for _ in 0:5]
  for i in 1:n, j in 1:n, (beta, value) in matrix[i, j].terms
    degree = sum(beta .* chern_coordinates)
    @assert degree == 1 + codimensions[j] - codimensions[i]
    @assert 0 <= degree <= 5
    push!(powers[degree + 1], (i, j, value))
  end
  initial = zeros(Q, n)
  initial[top] = 1
  states = [initial]
  ordinary = Q[1]
  for degree in 1:order
    rhs = zeros(Q, n)
    for k in 1:min(5, degree)
      rhs .+= apply_coefficient(powers[k + 1], states[degree - k + 1])
    end
    term = rhs ./ degree
    state = copy(term)
    for _ in 1:4
      term = apply_coefficient(powers[1], term) ./ degree
      state .+= term
    end
    push!(states, state)
    push!(ordinary, state[top])
  end
  regularized = [value * factorial(big(d - 1)) for (d, value) in enumerate(ordinary)]
  (; ordinary=rational_text.(ordinary), regularized=rational_text.(regularized))
end

"""
Run one complete matrix after warming the classical term and one nonzero curve.
Use a fresh process per engine and case. The single full sample includes any
remaining compilation; its measured compilation time is recorded separately.
The input TOML supplies the original full beta list and archived source hash.
"""
function run_case(backend, input_path, output)
  backend in ("old", "new") || error("backend must be old or new")
  input = TOML.parsefile(input_path)
  source = input["old_source"]
  source_hash = bytes2hex(sha256(read(source)))
  source_hash == input["old_source_sha256"] || error("archived engine hash mismatch")
  if backend == "old"
    Base.include(OldEngine, source)
  end
  engine = backend == "old" ? OldEngine : GromovWittenFlags
  # Including the historical source defines methods in a newer world. Dispatch
  # once across this boundary; the numerical loops retain normal specialization.
  Base.invokelatest(measure_case, engine, backend, input, output, source_hash)
end

"""Construct one target, measure its full matrix, and save exact reproducible outputs."""
function measure_case(engine, backend, input, output, source_hash)
  BLAS.set_num_threads(1)
  PFV.configure_caches!(; budget=128 * 1024^2)
  Semisimple.configure_caches!(; budget=128 * 1024^2)
  label = input["label"]
  progress(label, backend, "constructing")
  construction = @timed begin
    Z = PFV.zero_locus(label)
    X = PFV.ambient_variety(Z)
    @assert PFV.dimension(Z) == 4
    @assert PFV.zerolocus62_label(PFV.defining_bundle(Z)) == label
    _, betas = engine.fano_curve_classes(Z)
    @assert betas == input["betas"]
    data = engine._gw_data(X)
    dim, c_full, _ = engine._complete_intersection_data(data, PFV.defining_bundle(Z))
    @assert dim == 4
    c = Int[c_full[i] for i in PFV.marked_nodes(X)]
    @assert c == input["anticanonical_coordinates"]
  end
  warmup_betas = betas[1:min(2, length(betas))]
  progress(label, backend, "warming classical term and first nonzero curve")
  warmup = @timed engine.small_quantum_multiplication(Z; betas=warmup_betas)
  checkpoint(
    output * ".progress.toml",
    Dict(
      "label" => label, "backend" => backend, "phase" => "full quantum matrix",
      "warmup_seconds" => warmup.time, "warmup_betas" => warmup_betas,
    ),
  )
  GC.gc()
  progress(label, backend, "full quantum matrix; warmup $(warmup.time) s")
  timing = @timed engine.small_quantum_multiplication(Z; betas)
  result = timing.value
  progress(label, backend, "matrix complete in $(timing.time) s; period recurrence")
  codimensions = length.(result.basis)
  period_run = @timed period_coefficients(result.matrix, codimensions, c)
  period = period_run.value
  value = Dict(
    "label" => label, "backend" => backend, "julia_version" => string(VERSION),
    "new_source_commit" => input["new_source_commit"], "old_source_sha256" => source_hash,
    "partial_flag_varieties_version" => string(pkgversion(PFV)),
    "semisimple_version" => string(pkgversion(Semisimple)),
    "threads" => Threads.nthreads(), "gc_threads" => Threads.ngcthreads(),
    "blas_threads" => BLAS.get_num_threads(), "cpu" => Sys.CPU_NAME,
    "construction_seconds" => construction.time,
    "warmup_seconds" => warmup.time, "warmup_betas" => warmup_betas,
    "matrix_seconds" => timing.time, "matrix_allocated_bytes" => timing.bytes,
    "matrix_gc_seconds" => timing.gctime, "matrix_compile_seconds" => timing.compile_time,
    "period_seconds" => period_run.time, "observed_utc" => string(now(UTC)),
    "betas" => betas, "anticanonical_coordinates" => c,
    "matrix" => [[polynomial_text(p) for p in row] for row in eachrow(result.matrix)],
    "metric" => [[rational_text(x) for x in row] for row in eachrow(result.metric)],
    "basis_codimensions" => codimensions,
    "basis_keys" => [basis_key(data, w, input["ambient_factors"]) for w in result.basis],
    "grading" => rational_text.(diag(result.grading)),
    "unregularized_coefficients" => period.ordinary,
    "regularized_coefficients" => period.regularized,
  )
  checkpoint(output, value)
  progress(label, backend, "saved $output")
  value
end

end

if abspath(PROGRAM_FILE) == @__FILE__
  length(ARGS) == 3 || error("usage: hpc_comparison.jl old|new INPUT.toml OUTPUT.toml")
  HPCComparison.run_case(ARGS...)
end
