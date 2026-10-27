#!/usr/bin/env julia
# Keep the formatter version here as the single source of truth for local/CI use.
const JULIAFORMATTER_VERSION = "2.9.3"

using Pkg
Pkg.activate(; temp=true)
Pkg.add(
  Pkg.PackageSpec(; name="JuliaFormatter", version=JULIAFORMATTER_VERSION); io=devnull
)
using JuliaFormatter

const PACKAGE_ROOT = dirname(@__DIR__)
files = String[]
for directory in ("src", "test", "examples", "benchmark", "devtools", "docs")
  for (root, dirs, names) in walkdir(joinpath(PACKAGE_ROOT, directory))
    # Local benchmark archives can contain immutable historical Julia sources.
    filter!(directory -> directory ∉ ("build", "hpc-local"), dirs)
    append!(files, [joinpath(root, name) for name in names if endswith(name, ".jl")])
  end
end
check = "--check" in ARGS
results = [format(file; overwrite=(!check)) for file in files]
if check && !all(results)
  foreach(println, files[.!results])
  error("Formatting check failed; run julia devtools/format.jl")
end
println(check ? "Formatting verified." : "Formatting complete.")
