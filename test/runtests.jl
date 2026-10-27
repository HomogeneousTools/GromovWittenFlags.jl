using Test
using GromovWittenFlags
using PartialFlagVarieties
using LinearAlgebra: Diagonal, I
using Combinatorics: permutations, combinations
using Aqua

const GWF = GromovWittenFlags
const QQ = Rational{BigInt}

@testset "GromovWittenFlags.jl" begin
  @testset "Package quality" begin
    Aqua.test_all(GWF; ambiguities=(; recursive=false))
    missing = Symbol[]
    for name in names(GWF; all=true, imported=false)
      startswith(string(name), "_") || continue
      value = getfield(GWF, name)
      owned = (value isa Function || value isa Type) && parentmodule(value) === GWF
      (owned || isconst(GWF, name)) || continue
      # `doc` renders REPL help; `hasdoc` checks metadata without loading REPL.
      Base.Docs.hasdoc(GWF, name) || push!(missing, name)
    end
    @test isempty(missing)
  end
  include("schubert.jl")
  include("novikov.jl")
  include("trees.jl")
  include("invariants.jl")
  include("localization.jl")
  include("quantum.jl")
  include("paper_examples.jl")
  if get(ENV, "GWFLAGS_EXTENDED_TESTS", "false") == "true"
    include("kuechle_c5_quantum.jl")
    include("kuechle_b8_quantum.jl")
  end
end
