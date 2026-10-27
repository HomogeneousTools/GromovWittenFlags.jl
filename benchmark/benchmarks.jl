# Run from any checkout with `julia --project=. benchmark/benchmarks.jl`.
# Each workload is warmed once; reported medians exclude compilation.
"""Developer workloads and serial measurement helpers; inclusion does not start a run."""
module GromovWittenBenchmarks

using GromovWittenFlags
using PartialFlagVarieties

"""
Warm a callable once, then print median time and allocated bytes from odd samples.
Check each result against the warm reference outside timing. Inclusion defines
the driver without executing it; the CLI guard below starts a complete run.
"""
function measure(label, f::F; samples=3) where {F}
  samples > 0 && isodd(samples) || throw(ArgumentError("use a positive odd sample count"))
  reference = f()
  times = Float64[]
  bytes = Int[]
  for _ in 1:samples
    GC.gc()
    run = @timed f()
    run.value == reference || error("benchmark result changed: $label")
    push!(times, run.time)
    push!(bytes, run.bytes)
  end
  sort!(times)
  sort!(bytes)
  middle = samples ÷ 2 + 1
  println(label, '\t', times[middle], '\t', bytes[middle])
  flush(stdout)
  reference
end

"""Exercise the cover-degree worker repeatedly without retaining its borrowed buffers."""
function degree_search_workload(data; repetitions=20_000)
  count = Ref(0)
  count_assignment = _ -> (count[] += 1)
  roots, beta = [1, 1, 1], [5]
  for _ in 1:repetitions
    GromovWittenFlags._foreach_degree_assignment(count_assignment, data, roots, beta)
  end
  count[]
end

"""Exercise finite curve-vector enumeration, consuming each independently owned result."""
function curve_search_workload(; repetitions=10_000)
  count = 0
  for _ in 1:repetitions
    count += length(GromovWittenFlags._compact_beta_vectors([1, 2, 3], 8))
  end
  count
end

"""Run exact invariant, quantum-matrix, and enumeration workloads in a fixed order."""
function benchmarks()
  samples = parse(Int, get(ENV, "GWFLAGS_BENCH_SAMPLES", "3"))
  samples > 0 && isodd(samples) || error("use a positive odd sample count")
  println("workload\tmedian_seconds\tmedian_bytes")
  plane = projective_space(2)
  pt = point_class(plane)
  measure("P2 conics", () -> gromov_witten_invariant(plane, fill(pt, 5), [2]); samples)
  grassmannian = Gr(3, 7)
  data = GromovWittenFlags._gw_data(grassmannian)
  measure("degree search", () -> degree_search_workload(data); samples)
  measure("curve search", curve_search_workload; samples)
  E =
    exterior_power(dual(universal_subbundle(grassmannian)), 2) +
    exterior_power(universal_quotient_bundle(grassmannian), 3) +
    line_bundle(grassmannian, 1)
  c5 = zero_locus(E)
  measure("c5 classical", () -> small_quantum_multiplication(c5; betas=[[0]]); samples)
  measure(
    "c5 degree 2", () -> small_quantum_multiplication(c5; betas=[[0], [1], [2]]); samples
  )
  measure(
    "c5 degree 3", () -> small_quantum_multiplication(c5; betas=[[0], [1], [2], [3]]);
    samples,
  )
  quartic_ambient = projective_space(5)
  measure(
    "quartic fourfold",
    () -> small_quantum_multiplication(zero_locus(line_bundle(quartic_ambient, 4)));
    samples,
  )
  product_ambient = product(projective_space(2), projective_space(4))
  k329 = zero_locus(
    line_bundle(product_ambient, [1, 2]) + line_bundle(product_ambient, [1, 1])
  )
  measure("K3-29", () -> small_quantum_multiplication(k329); samples)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
  GromovWittenBenchmarks.benchmarks()
end
