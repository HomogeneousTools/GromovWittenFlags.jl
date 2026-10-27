include("GKMtoolsBenchmarks.jl")

if abspath(PROGRAM_FILE) == @__FILE__
  GKMtoolsBenchmarks.TOML.print(stdout, GKMtoolsBenchmarks.weyl_checks())
end
