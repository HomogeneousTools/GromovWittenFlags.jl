using Test
using TOML
include("GKMtoolsBenchmarks.jl")

const Bench = GKMtoolsBenchmarks

@testset "Comparison configuration and checkpoints" begin
  @test_throws ArgumentError Bench.benchmark_config(["misspelled"], Dict{String,String}())
  @test_throws ArgumentError Bench.benchmark_config(["cubic"],
    Dict("GWFLAGS_GKM_BASIS" => "misspelled"))
  @test_throws ArgumentError Bench.benchmark_config(["cubic"],
    Dict("GWFLAGS_BACKENDS" => "misspelled"))
  @test_throws ArgumentError Bench.backend_function("misspelled")
  @test !Bench.benchmark_config(
    ["cubic"],
    Dict(
      "GWFLAGS_GKM_PREPARED" => "true", "GWFLAGS_BACKENDS" => "GromovWittenFlags"),
  ).use_prepared

  mktempdir() do directory
    path = joinpath(directory, "checkpoint.toml")
    expected = Dict("completed" => 1)
    Bench.save_results(path, expected)
    @test TOML.parsefile(path) == expected
    @test_throws ErrorException Bench.save_results(path, Dict("unsupported" => identity))
    @test TOML.parsefile(path) == expected
    @test readdir(directory) == ["checkpoint.toml"]
  end
end

@testset "Weyl word transport" begin
  checks = Bench.weyl_checks()
  @test sum(row["simple_root_action_checks"] for row in checks["cases"]) == 256
end

@testset "Native matrix and prepared adapters" begin
  mktempdir() do directory
    config = Bench.benchmark_config(
      ["cubic", "K3-29"],
      Dict(
        "GWFLAGS_BENCH_SAMPLES" => "1", "GWFLAGS_GKM_PREPARED" => "true",
        "GWFLAGS_TOTAL_DEGREE" => "1",
        "GWFLAGS_BENCH_OUTPUT" => joinpath(directory, "results.toml")),
    )
    # The driver compares every coefficient with the independent exact result.
    # K3-29 exercises a nonconstant anticanonical degree vector, unlike Picard rank one.
    results = Bench.run_fourfolds(config)
    @test TOML.parsefile(config.output) == results
    @test length(results["cases"]) == 2
    @test all(row -> length(row["measurements"]) == 3, results["cases"])
  end
end
