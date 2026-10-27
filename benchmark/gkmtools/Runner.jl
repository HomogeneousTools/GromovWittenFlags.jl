"""Native model names accepted by the comparison driver, in its default run order."""
const FOURFOLD_NAMES = ("cubic", "Gushel-Mukai", "quartic", "quadric-cubic",
  "three-quadrics", "K3-29", "c5", "b8")

"""Public matrix APIs and the explicitly prepared integration adapter."""
const BACKEND_NAMES = ("GromovWittenFlags", "GKMtools", "GKMtools-prepared")

"""
    benchmark_config(args=ARGS, environment=ENV)

Parse and validate CLI model names and environment options before constructing
geometry. Return an owned named tuple. Accepting the environment as an argument
makes configurations testable without changing the process's global `ENV`.
Explicit `GWFLAGS_BACKENDS` takes precedence over the prepared-backend shortcut.
"""
function benchmark_config(args=ARGS, environment=ENV)
  samples = parse(Int, get(environment, "GWFLAGS_BENCH_SAMPLES", "3"))
  samples > 0 && isodd(samples) || throw(ArgumentError("use a positive odd sample count"))
  maximum_degree = parse(Int, get(environment, "GWFLAGS_KUECHLE_DEGREE", "3"))
  maximum_degree >= 0 || throw(ArgumentError("the degree cap must be nonnegative"))
  basis_lift = Symbol(get(environment, "GWFLAGS_GKM_BASIS", "billey"))
  basis_lift in (:billey, :schur) ||
    throw(ArgumentError("the basis lift must be billey or schur"))
  total_degree_text = get(environment, "GWFLAGS_TOTAL_DEGREE", "")
  total_degree = isempty(total_degree_text) ? nothing : parse(Int, total_degree_text)
  isnothing(total_degree) || total_degree >= 0 ||
    throw(ArgumentError("the degree cap must be nonnegative"))
  split_degrees = parse(Bool, get(environment, "GWFLAGS_SPLIT_DEGREES", "false"))
  add_prepared = parse(Bool, get(environment, "GWFLAGS_GKM_PREPARED", "false"))
  default_backends = add_prepared ? join(BACKEND_NAMES, ',') : join(BACKEND_NAMES[1:2], ',')
  backends = String.(
    strip.(split(get(environment, "GWFLAGS_BACKENDS", default_backends), ','))
  )
  all(in(BACKEND_NAMES), backends) || throw(ArgumentError("unrecognized backend"))
  use_prepared = "GKMtools-prepared" in backends
  names = isempty(args) ? collect(FOURFOLD_NAMES) : String.(args)
  all(in(FOURFOLD_NAMES), names) || throw(ArgumentError("unrecognized fourfold"))
  output = get(
    environment, "GWFLAGS_BENCH_OUTPUT", joinpath(@__DIR__, "results-local.toml")
  )
  (; samples, maximum_degree, total_degree, basis_lift, split_degrees,
    use_prepared, backends, names, output)
end

"""Resolve a validated CLI label to a computation; unknown labels never select a fallback."""
function backend_function(backend)
  backend == "GromovWittenFlags" && return flag_operator
  backend == "GKMtools" && return gkm_operator
  backend == "GKMtools-prepared" && return gkm_prepared_operator
  throw(ArgumentError("unrecognized backend: $backend"))
end

"""Choose native output conversion and specialize the measurement worker on both callables."""
function measure_backend(context, backend, expected; samples)
  compute = backend_function(backend)
  extract = if backend == "GromovWittenFlags"
    result -> flag_coefficients(context, result)
  else
    gkm_coefficients
  end
  measure_computation(compute, extract, context, expected; samples, backend)
end

"""
    measure_computation(compute, extract, context, expected; samples, backend)

Warm once, then record serial elapsed times and Julia-reported allocated bytes.
Validate every complete coefficient list outside timing. `compute` and `extract`
are specialized callables; string labels and serialization dictionaries stay out
of the timed kernel. The returned dictionary is a TOML reporting boundary.
"""
function measure_computation(compute::F, extract::E, context, expected;
  samples, backend) where {F,E}
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

"""
    save_results(path, results)

Serialize a completed checkpoint to a temporary file in the output directory,
then move it into place. Serialization errors leave the preceding checkpoint
intact. Temporary files are cleaned up on exceptions; the caller owns `results`.
"""
function save_results(path, results)
  mktemp(dirname(path)) do temporary, io
    TOML.print(io, results)
    close(io)
    mv(temporary, path; force=true)
  end
end

"""Retain corresponding native curve classes when selecting individual degrees."""
function select_degrees(context, indices)
  merge(context, (; betas=context.betas[indices], gkm_betas=context.gkm_betas[indices]))
end

"""Measure one backend, append its completed record, and persist failures before rethrowing."""
function measure_case!(row, context, backend, expected, results, config)
  println("Measuring ", context.name, " ", context.betas, " / ", backend)
  flush(stdout)
  started = time_ns()
  measurement = try
    measure_backend(context, backend, expected; samples=config.samples)
  catch exception
    push!(
      row["measurements"],
      Dict("backend" => backend,
        "status" => exception isa InterruptException ? "interrupted" : "error",
        "elapsed_seconds" => (time_ns() - started) / 1e9,
        "exception" => sprint(showerror, exception)),
    )
    save_results(config.output, results)
    rethrow()
  end
  push!(row["measurements"], measurement)
  save_results(config.output, results)
  println(backend, '\t', measurement["median_seconds"], '\t', measurement["median_bytes"])
  flush(stdout)
  nothing
end

"""Prepare one geometry, then append whole-operator or per-degree measurements to `results`."""
function benchmark_fourfold!(results, name, config)
  println("Preparing ", name)
  flush(stdout)
  setup = @timed prepare_fourfold(name;
    maximum_degree=config.maximum_degree, basis_lift=config.basis_lift)
  context, preparation_seconds = if config.use_prepared
    preparation = @timed prepare_gkm_operator(setup.value)
    merge(setup.value, (; prepared=preparation.value)), preparation.time
  else
    setup.value, 0.0
  end
  selected = if isnothing(config.total_degree)
    context
  else
    select_degrees(context,
      findall(beta -> sum(beta) <= config.total_degree, context.betas))
  end
  groups = if config.split_degrees
    [select_degrees(selected, [i]) for i in eachindex(selected.betas)]
  else
    [selected]
  end
  for group in groups
    expected = flag_coefficients(group, flag_operator(group))
    row = Dict{String,Any}("name" => name, "setup_seconds" => setup.time,
      "betas" => group.betas, "basis_size" => length(group.basis),
      "measurements" => Dict{String,Any}[])
    config.use_prepared && (row["gkm_preparation_seconds"] = preparation_seconds)
    push!(results["cases"], row)
    for backend in config.backends
      measure_case!(row, group, backend, expected, results, config)
    end
  end
  nothing
end

"""Run the selected comparison and return its saved, serializable measurement records."""
function run_fourfolds(config=benchmark_config())
  results = Dict{String,Any}("julia" => string(VERSION), "threads" => Threads.nthreads(),
    "timestamp_utc" => string(Dates.now(Dates.UTC)),
    "compiled_modules" => Base.JLOptions().use_compiled_modules,
    "packages" => Dict(
      string(nameof(mod)) => string(Base.pkgversion(mod))
      for mod in (GW, PF, SS, GK, OS)
    ),
    "gkmtools_revision" => "a19904bdb5b9c0158163b356dffb1cb79fd7a3c8",
    "samples" => config.samples, "basis_lift" => string(config.basis_lift),
    "split_degrees" => config.split_degrees, "prepared_backend" => config.use_prepared,
    "maximum_degree" => config.maximum_degree, "backends" => config.backends,
    "cases" => Dict{String,Any}[])
  isnothing(config.total_degree) || (results["total_degree"] = config.total_degree)
  for name in config.names
    benchmark_fourfold!(results, name, config)
  end
  results
end
