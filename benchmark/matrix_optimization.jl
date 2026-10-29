"""Local matrix benchmarks checked against twenty exact archived references."""
module MatrixOptimization

using GromovWittenFlags
using PartialFlagVarieties
using LinearAlgebra: BLAS, diag
using TOML
import Semisimple

include("hpc_comparison.jl")

"""Active uptime excludes macOS system suspension from an interrupted measurement."""
function active_seconds()
  Sys.isapple() || return time_ns() / 1e9
  info = Ref{NTuple{2,UInt32}}((0, 0))
  ccall(:mach_timebase_info, Cint, (Ref{NTuple{2,UInt32}},), info) == 0 ||
    error("mach_timebase_info failed")
  ticks = ccall(:mach_absolute_time, UInt64, ())
  numerator, denominator = info[]
  Float64(ticks) * (Float64(numerator) / Float64(denominator)) / 1e9
end

"""Canonical full matrix, metric, basis, grading, curves, and degree-100 periods."""
function canonical_result(result, data, record)
  codimensions = length.(result.basis)
  reference = record["reference"]
  period = HPCComparison.period_coefficients(result.matrix, codimensions,
    reference["anticanonical_coordinates"], 100)
  Dict(
    "matrix" => [[HPCComparison.polynomial_text(p) for p in row]
     for row in eachrow(result.matrix)],
    "metric" => [[HPCComparison.rational_text(x) for x in row]
     for row in eachrow(result.metric)],
    "grading" => HPCComparison.rational_text.(diag(result.grading)),
    "basis_keys" => [HPCComparison.basis_key(data, w, record["ambient_factors"])
     for w in result.basis],
    "basis_codimensions" => codimensions,
    "betas" => reference["betas"],
    "anticanonical_coordinates" => reference["anticanonical_coordinates"],
    "unregularized_coefficients" => period.ordinary,
    "regularized_coefficients" => period.regularized,
  )
end

"""Time complete matrix calls, then compare every exact output outside the timer."""
function run_case(record; warm=false)
  label = record["label"]
  Z = zero_locus(label)
  X = ambient_variety(Z)
  data = GromovWittenFlags._gw_data(X)
  dimension(Z) == 4 || error("expected a fourfold: $label")
  dimension(X) == record["ambient_dimension"] || error("ambient dimension changed")
  _, betas = fano_curve_classes(Z)
  reference = record["reference"]
  betas == reference["betas"] || error("curve classes changed: $label")
  _, c_full, _ = GromovWittenFlags._complete_intersection_data(data, defining_bundle(Z))
  c = [c_full[node] for node in marked_nodes(X)]
  c == reference["anticanonical_coordinates"] || error("first Chern class changed")
  rows = Any[]
  for sample in 1:(warm ? 2 : 1)
    GC.gc()
    start = active_seconds()
    timing = @timed small_quantum_multiplication(Z; betas)
    elapsed = active_seconds() - start
    canonical_result(timing.value, data, record) == reference ||
      error("exact reference mismatch: $label")
    push!(
      rows,
      Dict("sample" => sample, "active_seconds" => elapsed,
        "seconds_including_suspension" => timing.time,
        "compile_seconds" => timing.compile_time, "julia_allocation_bytes" => timing.bytes,
        "exact_match" => true),
    )
    println(label, " sample ", sample, ": ", round(elapsed; digits=3), " s, exact match")
    flush(stdout)
  end
  Dict("label" => label, "samples" => rows)
end

"""Run all archived cases or the supplied labels; `--warm` adds a second matrix call."""
function main(args=ARGS)
  warm = "--warm" in args
  labels = filter(!=("--warm"), args)
  archive = TOML.parsefile(joinpath(@__DIR__, "matrix-optimization", "references.toml"))
  cases = archive["cases"]
  all(label -> any(case -> case["label"] == label, cases), labels) ||
    error("unknown model label")
  BLAS.set_num_threads(1)
  PartialFlagVarieties.configure_caches!(; budget=128 * 1024^2)
  Semisimple.configure_caches!(; budget=128 * 1024^2)
  selected = filter(case -> isempty(labels) || case["label"] in labels, cases)
  output = joinpath(@__DIR__, "matrix-optimization", "local-results.toml")
  result = Dict("julia_version" => string(VERSION), "cpu" => Sys.CPU_NAME,
    "threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
    "cases" => Any[])
  for case in selected
    push!(result["cases"], run_case(case; warm))
    HPCComparison.checkpoint(output, result)
  end
  println("Saved ", output)
  result
end

end

if abspath(PROGRAM_FILE) == @__FILE__
  MatrixOptimization.main()
end
