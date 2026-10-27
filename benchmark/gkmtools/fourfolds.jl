using Random
using TOML
using Dates
include("Models.jl")
include("PreparedOperator.jl")

"""Warm once, then record serial samples; validate every output outside the timed section."""
function measure_backend(context, backend, expected; samples)
  compute = if backend == "GromovWittenFlags"
    flag_operator
  elseif backend == "GKMtools"
    gkm_operator
  else
    gkm_prepared_operator
  end
  extract = if backend == "GromovWittenFlags"
    result -> flag_coefficients(context, result)
  else
    gkm_coefficients
  end
  Random.seed!(20261027)
  warm = @timed compute(context)
  extract(warm.value) == expected ||
    error("$backend coefficient mismatch on $(context.name)")
  elapsed, allocated = Float64[], Int[]
  for sample in 1:samples
    Random.seed!(20261027 + sample)
    GC.gc()
    result = @timed compute(context)
    extract(result.value) == expected ||
      error("$backend coefficient mismatch on $(context.name)")
    push!(elapsed, result.time)
    push!(allocated, result.bytes)
  end
  middle = samples ÷ 2 + 1
  Dict("backend" => backend, "warm_seconds" => warm.time,
    "seconds" => elapsed, "bytes" => allocated,
    "median_seconds" => sort(elapsed)[middle], "median_bytes" => sort(allocated)[middle])
end

"""Persist completed rows immediately, so a costly case cannot erase earlier measurements."""
function save_results(path, results)
  open(path, "w") do io
    TOML.print(io, results)
  end
end

"""Retain corresponding native curve classes when selecting individual degrees."""
function select_degrees(context, indices)
  merge(context, (; betas=context.betas[indices], gkm_betas=context.gkm_betas[indices]))
end

"""Benchmark common coefficient matrices, optionally one curve class at a time."""
function main()
  samples = parse(Int, get(ENV, "GWFLAGS_BENCH_SAMPLES", "3"))
  samples > 0 && isodd(samples) || error("use a positive odd sample count")
  maximum_degree = parse(Int, get(ENV, "GWFLAGS_KUECHLE_DEGREE", "3"))
  maximum_degree >= 0 || error("the degree cap must be nonnegative")
  basis_lift = Symbol(get(ENV, "GWFLAGS_GKM_BASIS", "billey"))
  total_degree_text = get(ENV, "GWFLAGS_TOTAL_DEGREE", "")
  total_degree = isempty(total_degree_text) ? nothing : parse(Int, total_degree_text)
  isnothing(total_degree) || total_degree >= 0 ||
    error("the degree cap must be nonnegative")
  split_degrees = parse(Bool, get(ENV, "GWFLAGS_SPLIT_DEGREES", "false"))
  use_prepared = parse(Bool, get(ENV, "GWFLAGS_GKM_PREPARED", "false"))
  default_backends = if use_prepared
    "GromovWittenFlags,GKMtools,GKMtools-prepared"
  else
    "GromovWittenFlags,GKMtools"
  end
  backends = strip.(split(get(ENV, "GWFLAGS_BACKENDS", default_backends), ','))
  all(
    backend -> backend in ("GromovWittenFlags", "GKMtools", "GKMtools-prepared"), backends
  ) ||
    error("unrecognized backend")
  use_prepared |= "GKMtools-prepared" in backends
  names = if isempty(ARGS)
    ["cubic", "Gushel-Mukai", "quartic", "quadric-cubic",
      "three-quadrics", "K3-29", "c5", "b8"]
  else
    ARGS
  end
  output = get(ENV, "GWFLAGS_BENCH_OUTPUT", joinpath(@__DIR__, "results-local.toml"))
  results = Dict{String,Any}("julia" => string(VERSION), "threads" => Threads.nthreads(),
    "timestamp_utc" => string(Dates.now(Dates.UTC)),
    "compiled_modules" => Base.JLOptions().use_compiled_modules,
    "packages" => Dict(
      string(nameof(mod)) => string(Base.pkgversion(mod))
      for mod in (GW, PF, SS, GK, OS)
    ),
    "gkmtools_revision" => "a19904bdb5b9c0158163b356dffb1cb79fd7a3c8",
    "samples" => samples, "cases" => Any[])
  results["basis_lift"] = string(basis_lift)
  results["split_degrees"] = split_degrees
  results["prepared_backend"] = use_prepared
  for name in names
    println("Preparing ", name)
    flush(stdout)
    setup = @timed prepare_fourfold(name; maximum_degree, basis_lift)
    context = setup.value
    preparation_seconds = 0.0
    if use_prepared
      preparation = @timed prepare_gkm_operator(context)
      context = merge(context, (; prepared=preparation.value))
      preparation_seconds = preparation.time
    end
    if !isnothing(total_degree)
      context = select_degrees(
        context, findall(beta -> sum(beta) <= total_degree, context.betas)
      )
    end
    groups = if split_degrees
      [select_degrees(context, [i]) for i in eachindex(context.betas)]
    else
      [context]
    end
    for context in groups
      expected = flag_coefficients(context, flag_operator(context))
      row = Dict{String,Any}("name" => name, "setup_seconds" => setup.time,
        "betas" => context.betas, "basis_size" => length(context.basis),
        "measurements" => Any[])
      use_prepared && (row["gkm_preparation_seconds"] = preparation_seconds)
      push!(results["cases"], row)
      for backend in backends
        println("Measuring ", name, " ", context.betas, " / ", backend)
        flush(stdout)
        started = time_ns()
        measurement = try
          measure_backend(context, backend, expected; samples)
        catch exception
          push!(
            row["measurements"],
            Dict("backend" => backend,
              "status" => exception isa InterruptException ? "interrupted" : "error",
              "elapsed_seconds" => (time_ns() - started) / 1e9,
              "exception" => sprint(showerror, exception)),
          )
          save_results(output, results)
          rethrow()
        end
        push!(row["measurements"], measurement)
        save_results(output, results)
        println(
          backend, '\t', measurement["median_seconds"], '\t', measurement["median_bytes"]
        )
        flush(stdout)
      end
    end
  end
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
