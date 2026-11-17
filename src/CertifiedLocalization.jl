"""Validate a prime once, before creating the specialized modular scalar type."""
function _cert_prime_type(p::Integer)
  2 <= p <= 2147483647 && Nemo.is_prime(p) ||
    throw(ArgumentError("certified moduli must be primes at most 2^31-1"))
  _CertPrime{Int(p)}
end

"""Hash primitive Julia inputs without relying on randomized object hashes."""
function _cert_digest(value)
  io = IOBuffer()
  Serialization.serialize(io, value)
  bytes2hex(SHA.sha256(take!(io)))
end

"""Content identity of the package sources and relevant dependency versions."""
function _cert_source_identity()
  files = sort(filter(f -> endswith(f, ".jl"), readdir(@__DIR__)))
  sources = [(f, bytes2hex(SHA.sha256(read(joinpath(@__DIR__, f))))) for f in files]
  push!(
    sources,
    ("Project.toml", bytes2hex(SHA.sha256(read(joinpath(@__DIR__, "..", "Project.toml"))))),
  )
  _cert_digest((string(VERSION), sources, string(Base.pkgversion(Nemo)),
    string(Base.pkgversion(parentmodule(PartialFlagVariety))),
    string(Base.pkgversion(parentmodule(WeightLatticeElem)))))
end

"""Serialize a checkpoint canonically; its digest detects truncated or altered records."""
function _cert_record_text(record)
  io = IOBuffer()
  TOML.print(io, record; sorted=true)
  String(take!(io))
end

"""Atomically replace one complete checkpoint on the same filesystem."""
function _cert_store(path, record, identity)
  path === nothing && return record
  mkpath(dirname(path))
  stored = merge(record, Dict("identity" => identity, "format_version" => 1))
  stored["checksum"] = bytes2hex(SHA.sha256(_cert_record_text(stored)))
  temporary, io = mktemp(dirname(path))
  try
    write(io, _cert_record_text(stored))
    flush(io)
    close(io)
    mv(temporary, path; force=true)
  finally
    isopen(io) && close(io)
    isfile(temporary) && rm(temporary)
  end
  record
end

"""Reject incompatible or corrupt checkpoints rather than silently reusing them."""
function _cert_load(path, identity)
  (path === nothing || !isfile(path)) && return nothing
  record = TOML.parsefile(path)
  get(record, "identity", "") == identity && get(record, "format_version", 0) == 1 ||
    throw(
      ArgumentError("checkpoint identity mismatch: $path; use a new checkpoint directory")
    )
  checksum = pop!(record, "checksum", "")
  checksum == bytes2hex(SHA.sha256(_cert_record_text(record))) ||
    throw(ArgumentError("corrupt checkpoint: $path"))
  record
end

"""Join a checkpoint path while supporting computations without disk checkpoints."""
_cert_path(directory, name) = directory === nothing ? nothing : joinpath(directory, name)

"""One process-local exact context; no FLINT objects cross process boundaries."""
const _CERT_CONTEXT = Ref{Any}(nothing)

"""Rebuild torus geometry locally and reuse it for a batch of primes."""
function _cert_context(spec, beta, seed)
  key = (spec.identity, beta, seed)
  cached = _CERT_CONTEXT[]
  cached !== nothing && cached[1] == key && return cached[2]
  _CERT_CONTEXT[] = nothing
  GC.gc()
  data = _gw_data(spec.variety)
  DT = dynkin_type(spec.variety)
  weights =
    spec.weights === nothing ? nothing : [WeightLatticeElem(DT, w) for w in spec.weights]
  evaluation = _gw_evaluation(data, seed, QQFieldElem; weights)
  geometry = _branch_geometry(evaluation, beta)
  table = _restriction_table(evaluation, spec.classes)
  context = (; evaluation, geometry, table, pairs=spec.pairs, beta)
  _CERT_CONTEXT[] = (key, context)
  context
end

"""Compute all requested pairs with shared marked-source recursions."""
function _cert_pair_values(context, ::Type{T}) where {T}
  evaluation = _cert_convert_evaluation(context.evaluation, T)
  geometry = _cert_convert_geometry(context.geometry, T)
  table = T.(context.table)
  sources = _branch_pair_sources(context.pairs)
  assignments = [j in sources ? j : i for (i, j) in context.pairs]
  values = zeros(T, length(context.pairs))
  for source in sources
    roots = _branch_root_values(evaluation, table[:, source:source], context.beta; geometry)
    for (index, (i, j)) in enumerate(context.pairs)
      assignments[index] == source || continue
      other = i == source ? j : i
      values[index] = sum(table[:, other] .* roots)
    end
    GC.gc()
  end
  values
end

"""Execute one bound or prime job, returning only primitive serializable values."""
function _cert_job(job)
  started = time_ns()
  try
    context = _cert_context(job.spec, job.beta, job.seed)
    record = if job.kind == :bound
      denominator = _cert_denominator_certificate(context.evaluation, job.beta)
      D = parse(BigInt, denominator["denominator_bound"])
      bounds = _cert_pair_values(context, _CertBound)
      scaled = [value * _CertBound(D) for value in bounds]
      Dict{String,Any}("kind" => "bound", "denominator" => denominator,
        "upper" => string.(_cert_bound_integer_ceiling.(scaled)),
        "mantissas" => string.([v.mantissa for v in scaled]),
        "exponents" => [v.exponent for v in scaled])
    else
      values = _cert_pair_values(context, _cert_prime_type(job.prime))
      Dict{String,Any}("kind" => "prime", "prime" => job.prime,
        "accepted" => true, "residues" => [v.value for v in values])
    end
    record["seed"] = job.seed
    record["beta"] = collect(job.beta)
    record["seconds"] = (time_ns() - started) / 1e9
    record["worker_pid"] = getpid()
    GC.gc()
    record
  catch error
    if error isa _LocalizationPole && job.kind == :prime
      return Dict{String,Any}("kind" => "prime", "prime" => job.prime,
        "seed" => job.seed, "beta" => collect(job.beta), "accepted" => false,
        "reason" => "singular specialization in this finite field",
        "seconds" => (time_ns() - started) / 1e9, "worker_pid" => getpid())
    end
    rethrow()
  end
end

"""Own a bounded process pool, or use caller-owned Distributed worker IDs."""
function _cert_with_pool(f, workers; worker_heap_size=nothing)
  owned = Int[]
  ids = if workers isa Integer
    workers >= 0 || throw(ArgumentError("workers must be nonnegative"))
    if workers > 0
      project = Base.active_project()
      project === nothing &&
        throw(ArgumentError("parallel certification requires an active Julia project"))
      flags = ["--project=$(dirname(project))", "--threads=1", "--startup-file=no"]
      worker_heap_size === nothing || push!(flags, "--heap-size-hint=$worker_heap_size")
      append!(
        owned,
        Distributed.addprocs(workers; exeflags=Cmd(flags),
          env=["OPENBLAS_NUM_THREADS" => "1", "OMP_NUM_THREADS" => "1"]),
      )
    end
    owned
  else
    Int.(collect(workers))
  end
  try
    length(unique(ids)) == length(ids) && all(i -> i != Distributed.myid(), ids) ||
      throw(ArgumentError("worker IDs must be distinct remote processes"))
    for id in ids
      Distributed.remotecall_fetch(Core.eval, id, Main, :(using GromovWittenFlags))
    end
    f(isempty(ids) ? nothing : Distributed.WorkerPool(ids))
  finally
    isempty(owned) || Distributed.rmprocs(owned)
    _CERT_CONTEXT[] = nothing
    GC.gc()
  end
end

"""Map independent jobs; a callback persists each completed result immediately."""
function _cert_map_jobs(consume, pool, jobs)
  if pool === nothing
    for job in jobs
      consume(job, _cert_job(job))
    end
  else
    # One asynchronous feeder per worker; completed jobs are checkpointed even
    # when another worker subsequently fails. Only the coordinator writes files.
    queue = Channel{Any}(length(jobs))
    foreach(job -> put!(queue, job), jobs)
    close(queue)
    @sync for _ in 1:min(length(jobs), length(Distributed.workers(pool)))
      @async for job in queue
        result = Distributed.remotecall_fetch(_cert_job, pool, job)
        consume(job, result)
      end
    end
  end
  nothing
end

"""Combine scaled rational residues and apply the exact bounded-integer test."""
function _cert_recover(records, denominator, upper)
  residue = fill(big(0), length(upper))
  modulus = big(1)
  used = Set{Int}()
  for record in records
    get(record, "accepted", false) || continue
    prime = record["prime"]
    prime in used && throw(ArgumentError("duplicate prime in a CRT certificate"))
    _cert_prime_type(prime)
    iszero(rem(denominator, prime)) &&
      throw(ArgumentError("prime divides denominator bound"))
    push!(used, prime)
    length(record["residues"]) == length(upper) ||
      throw(DimensionMismatch("CRT residue vector"))
    inverse = invmod(mod(modulus, prime), prime)
    for i in eachindex(upper)
      r = record["residues"][i]
      0 <= r < prime || throw(ArgumentError("noncanonical modular residue"))
      correction = mod((denominator * r - residue[i]) * inverse, prime)
      residue[i] += modulus * correction
    end
    modulus *= prime
  end
  values = _BigRational[]
  for i in eachindex(upper)
    lower = cld(-upper[i] - residue[i], modulus)
    higher = fld(upper[i] - residue[i], modulus)
    lower == higher || return nothing
    push!(values, (residue[i] + lower * modulus)//denominator)
  end
  (; values, modulus, residue, prime_count=length(used))
end

"""Choose enough descending good primes for the proved numerator bound."""
function _cert_prime_batch(denominator, upper, used; minimum=1)
  result = Int[]
  modulus = big(1)
  threshold = 2 * maximum(upper)
  p = 2147483647
  while length(result) < minimum || modulus <= threshold
    if !(p in used) && !iszero(rem(denominator, p)) && Nemo.is_prime(p)
      push!(result, p)
      modulus *= p
    end
    p -= 2
    p > 1073741824 || error("exhausted supported large primes")
  end
  result
end

"""Certify one dimension-zero Schubert pair block with two independent torus recoveries."""
function _cert_pair_block(data, classes, pairs, beta, weights, pool;
  checkpoint=nothing, torus_seeds=(0, 3), progress=nothing, max_primes=256)
  max_primes isa Integer && max_primes > 0 ||
    throw(ArgumentError("max_primes must be positive"))
  seeds = Int.(collect(torus_seeds))
  length(seeds) == 2 && length(unique(seeds)) == 2 && all(s -> 0 <= s < 32, seeds) ||
    throw(ArgumentError("choose two distinct torus seeds in 0:31"))
  _vertex_bound(beta)
  all(iszero, beta) && throw(ArgumentError("certification requires positive curve degree"))
  source = _cert_source_identity()
  primitive_weights =
    weights === nothing ? nothing : [collect(Int, coefficients(w)) for w in weights]
  identity = _cert_digest((source, marked_dynkin_type(data.variety), primitive_weights,
    collect(classes), collect(pairs), collect(beta), seeds))
  spec = (; identity, variety=data.variety, weights=primitive_weights,
    classes=collect(classes), pairs=collect(pairs))
  directory = _cert_path(checkpoint, "beta-" * join(beta, "-"))
  manifest = Dict{String,Any}("source_identity" => source, "beta" => collect(beta),
    "classes" => collect(classes), "pairs" => [collect(p) for p in pairs], "seeds" => seeds,
  )
  previous = _cert_load(_cert_path(directory, "manifest.toml"), identity)
  previous === nothing &&
    _cert_store(_cert_path(directory, "manifest.toml"), manifest, identity)
  started = time_ns()
  boundpath = _cert_path(directory, "bound.toml")
  bound = _cert_load(boundpath, identity)
  bound_reused = bound !== nothing
  if bound === nothing
    job = (; kind=:bound, spec, beta=collect(beta), seed=first(seeds), prime=0)
    _cert_map_jobs(pool, [job]) do _, record
      bound = _cert_store(boundpath, record, identity)
    end
  end
  denominator = parse(BigInt, bound["denominator"]["denominator_bound"])
  upper = parse.(BigInt, bound["upper"])
  bound["kind"] == "bound" && bound["beta"] == collect(beta) &&
  bound["seed"] == first(seeds) ||
    error("bound checkpoint does not belong to this job")
  denominator > 0 && length(upper) == length(pairs) && all(>=(0), upper) ||
    error("invalid denominator or bound vector")
  length(bound["mantissas"]) == length(upper) == length(bound["exponents"]) ||
    throw(DimensionMismatch("bound checkpoint vectors"))
  for (u, m, e) in zip(upper, bound["mantissas"], bound["exponents"])
    exact = parse(BigInt, m) * (big(2)//1)^e
    u == cld(numerator(exact), Base.denominator(exact)) || error("invalid bound checkpoint")
  end
  progress === nothing || progress((phase=:certified_bound, beta=collect(beta),
    bits=maximum(ndigits(u; base=2) for u in upper), reused=bound_reused))
  records = Dict(seed => Dict{Int,Any}() for seed in seeds)
  computed = Ref(0)
  reused = Ref(0)
  used = Set{Int}()
  primes = _cert_prime_batch(denominator, upper, used)
  recoveries = Dict{Int,Any}()
  while true
    union!(used, primes)
    length(used) <= max_primes ||
      error("max_primes exhausted without a unique exact recovery")
    jobs = Any[]
    for seed in seeds, prime in primes
      prime > sum(beta) || throw(ArgumentError("prime must exceed the total degree"))
      path = _cert_path(directory, "seed-$seed-prime-$prime.toml")
      record = _cert_load(path, identity)
      if record === nothing
        push!(jobs, (; kind=:prime, spec, beta=collect(beta), seed, prime))
      else
        record["kind"] == "prime" && record["beta"] == collect(beta) &&
        record["seed"] == seed && record["prime"] == prime ||
          error("prime checkpoint does not belong to this job")
        records[seed][prime] = record
        reused[] += 1
      end
    end
    _cert_map_jobs(pool, jobs) do job, record
      path = _cert_path(directory, "seed-$(job.seed)-prime-$(job.prime).toml")
      records[job.seed][job.prime] = _cert_store(path, record, identity)
      computed[] += 1
      progress === nothing || progress((phase=:certified_prime, beta=collect(beta),
        seed=job.seed, prime=job.prime, accepted=record["accepted"],
        seconds=record["seconds"]))
    end
    empty!(recoveries)
    for seed in seeds
      ordered = [records[seed][p] for p in sort(collect(keys(records[seed])); rev=true)]
      recovery = _cert_recover(ordered, denominator, upper)
      recovery === nothing || (recoveries[seed] = recovery)
    end
    length(recoveries) == 2 && break
    primes = _cert_prime_batch(
      denominator, [big(0)], used; minimum=max(1, length(jobs) ÷ 8)
    )
  end
  first_recovery, second_recovery = (recoveries[s] for s in seeds)
  first_recovery.values == second_recovery.values ||
    error("exact torus recoveries disagree")
  certificate = Dict{String,Any}("status" => "certified", "beta" => collect(beta),
    "identity" => identity, "source_identity" => source,
    "denominator_bound" => string(denominator), "upper" => string.(upper),
    "values" => string.(first_recovery.values), "seeds" => seeds,
    "moduli" => [string(recoveries[s].modulus) for s in seeds],
    "prime_counts" => [recoveries[s].prime_count for s in seeds],
    "new_prime_jobs" => computed[], "reused_prime_jobs" => reused[],
    "bound_reused" => bound_reused, "seconds" => (time_ns() - started) / 1e9)
  _cert_store(_cert_path(directory, "certificate.toml"), certificate, identity)
  (; values=first_recovery.values, certificate)
end

"""
    certified_pair_invariants(X, classes, pairs, beta; bundle=nothing,
        workers=0, checkpoint=nothing, torus_seeds=(0,3), max_primes=256)

Compute exact dimension-zero two-point invariants of integral Schubert classes.
`pairs` indexes `classes`. The result contains exact rational `values` and an
arithmetic `certificate`. Formal graph-factor denominators and rigorous integer
magnitude bounds make CRT recovery deterministic; integrality is never assumed.

`workers` is a number of owned local Distributed processes, or existing worker
IDs. Zero runs serially. Each process owns its arithmetic and geometry caches.
Completed prime jobs and bounds are atomically checkpointed as checksummed TOML.
Source, dependency, model, insertion, degree and torus identities must agree on
resume. One coordinator should own a checkpoint directory. Failed jobs leave
completed checkpoints usable. The two distinct exact torus specializations
must be nonsingular; select other seeds if a torus preflight reports a pole.
"""
function certified_pair_invariants(X::PartialFlagVariety, classes, pairs, beta;
  bundle=nothing, workers=0, worker_heap_size=nothing, checkpoint=nothing,
  torus_seeds=(0, 3), max_primes=256, progress=nothing)
  data = _gw_data(X)
  dim, c, weights = _complete_intersection_data(data, bundle)
  normalized = _normalize_beta(data, beta)
  indices = [_class_index(data, s) for s in classes]
  required = dim - 1 + _curve_degree(data, c, normalized)
  all(p -> length(p) == 2 && all(i -> 1 <= i <= length(indices), p), pairs) ||
    throw(ArgumentError("pairs must index the supplied classes"))
  all(p -> sum(data.codimensions[indices[i]] for i in p) == required, pairs) ||
    throw(ArgumentError("insertions must have dimension-zero virtual degree"))
  isempty(pairs) && throw(ArgumentError("at least one insertion pair is required"))
  _cert_with_pool(workers; worker_heap_size) do pool
    _cert_pair_block(data, indices, Tuple{Int,Int}.(pairs), normalized, weights, pool;
      checkpoint, torus_seeds, progress, max_primes)
  end
end
