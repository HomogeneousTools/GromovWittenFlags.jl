# Run from any checkout with `julia --project=. benchmark/benchmarks.jl`.
# Each workload is warmed once; reported medians exclude compilation.
using GromovWittenFlags
using PartialFlagVarieties

function measure(label, f; samples=3)
  reference = f()
  times = Float64[]
  bytes = Int[]
  for _ in 1:samples
    GC.gc()
    run = @timed f()
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

function benchmarks()
  samples = parse(Int, get(ENV, "GWFLAGS_BENCH_SAMPLES", "3"))
  samples > 0 && isodd(samples) || error("use a positive odd sample count")
  println("workload\tmedian_seconds\tmedian_bytes")
  X = projective_space(2)
  pt = point_class(X)
  measure("P2 conics", () -> gromov_witten_invariant(X, fill(pt, 5), [2]); samples)
  G = Gr(3, 7)
  E =
    exterior_power(dual(universal_subbundle(G)), 2) +
    exterior_power(universal_quotient_bundle(G), 3) + line_bundle(G, 1)
  Z = zero_locus(E)
  measure("c5 classical", () -> small_quantum_multiplication(Z; betas=[[0]]); samples)
  measure(
    "c5 degree 2", () -> small_quantum_multiplication(Z; betas=[[0], [1], [2]]); samples
  )
  measure(
    "c5 degree 3", () -> small_quantum_multiplication(Z; betas=[[0], [1], [2], [3]]);
    samples,
  )
  X = projective_space(5)
  measure(
    "quartic fourfold",
    () -> small_quantum_multiplication(zero_locus(line_bundle(X, 4)));
    samples,
  )
  X = product(projective_space(2), projective_space(4))
  Z = zero_locus(line_bundle(X, [1, 2]) + line_bundle(X, [1, 1]))
  measure("K3-29", () -> small_quantum_multiplication(Z); samples)
end

benchmarks()
