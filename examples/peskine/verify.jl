# Post-computation comparison. This fixture is never an input to compute.jl.
using Test
import TOML, Nemo
include("graded_algebra.jl")

function exact_rational(s)
  parts = split(replace(s, "//"=>"/"), "/")
  parse(BigInt, parts[1]) // (length(parts)==1 ? big(1) : parse(BigInt, parts[2]))
end

length(ARGS)>=1 || error("usage: verify.jl RUN_DIRECTORY [fresh|resume]")
directory = ARGS[1]
report = TOML.parsefile(joinpath(directory, "peskine.toml"))
fixture = TOML.parsefile(
  joinpath(@__DIR__, "fixtures", "peskine-previous-certificate.toml")
)
@testset "Complete regenerated Peskine computation" begin
  @test report["status"]=="complete"
  @test length(report["coefficient_blocks"])==6
  @test length(report["two_point_invariants"])==12
  for d in 0:5
    block = only(filter(x -> x["degree"]==d, report["coefficient_blocks"]))
    reference = only(filter(x -> x["degree"]==d, fixture["coefficients"]))
    @test [exact_rational.(r) for r in block["matrix"]] ==
      [exact_rational.(r) for r in reference["matrix"]]
  end
  R, t = Nemo.polynomial_ring(Nemo.QQ, "t")
  @test report["characteristic_polynomial_at_q_1"]==string((t-228)*(t+14)^2*(t+15)^3)
  @test report["minimal_polynomial_at_q_1"]==report["characteristic_polynomial_at_q_1"]
  @test exact_rational(report["krylov_determinant_at_q_1"])==32076//25
  @test exact_rational(report["full_degree_one_trace"]) == -633
  @test exact_rational.(report["line_surface_moments"]) == [1179, 534, 1728, 2448, 876]
  @test exact_rational(report["inherited_scalar"]) == -14
  @test exact_rational(report["vanishing_scalar"]) == -15
  @test sort([b["ambient_dimension"] for b in report["spectral_blocks"]])==[1, 2, 3]
  @test sort([b["full_dimension"] for b in report["spectral_blocks"]])==[1, 24, 35]
  @test length(report["multiplication_operators_at_q_1"])==6
  @test length(report["certification"])==5
  @test export_graded_algebra(directory)["status"]=="passed"
  for c in report["certification"]
    @test c["status"]=="certified"
    if length(ARGS)>1 && ARGS[2]=="fresh"
      @test !c["bound_reused"]
      @test c["reused_prime_jobs"]==0
      @test c["new_prime_jobs"]>0
    elseif length(ARGS)>1 && ARGS[2]=="resume"
      @test c["bound_reused"]
      @test c["new_prime_jobs"]==0
    end
  end
end
