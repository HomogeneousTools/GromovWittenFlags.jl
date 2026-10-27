"""
Isolated comparison tools for native GKMtools and GromovWittenFlags models.

Loading this module defines the drivers without running workloads. Configuration
is parsed at the entry point; model construction and timed arithmetic remain in
separate functions. No names or methods are added to either library's module.
"""
module GKMtoolsBenchmarks

using Dates
using Random
using TOML

include("Models.jl")
include("PreparedOperator.jl")
include("Runner.jl")
include("WeylChecks.jl")

end
