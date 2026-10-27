include("GKMtoolsBenchmarks.jl")

if abspath(PROGRAM_FILE) == @__FILE__
  GKMtoolsBenchmarks.run_fourfolds()
end
