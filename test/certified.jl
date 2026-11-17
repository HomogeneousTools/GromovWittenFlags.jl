import Nemo, TOML, Distributed

@testset "Certified integer and modular arithmetic" begin
  Q = Rational{BigInt}
  boundvalue(x) = big(x.mantissa)*(big(2)//1)^x.exponent
  samples = [Q(0), Q(1), Q(-1), big(1)//3, big(-7)//15, big(123456789)//77,
    (big(1)<<200)//13, big(-3)//(big(1)<<240)]
  @test isbitstype(GWF._CertBound)
  @test sizeof(GWF._CertBound)==32
  for a in samples, b in samples
    x, y = GWF._CertBound(a), GWF._CertBound(b)
    @test boundvalue(x)>=abs(a)
    @test GWF._CertBound(Nemo.QQ(a))==x
    @test boundvalue(x+y)>=abs(a+b)
    @test boundvalue(x-y)>=abs(a-b)
    @test boundvalue(x*y)>=abs(a*b)
    iszero(b) || @test boundvalue(x/y)>=abs(a/b)
  end
  T = GWF._cert_prime_type(2147483647)
  for a in -4:4, b in -4:4
    @test (T(a)+T(b)).value == mod(a+b, 2147483647)
    @test (T(a)*T(b)).value == mod(a*b, 2147483647)
    iszero(b) || @test T(a)/T(b) == T(a//b)
  end
  @test_throws GWF._LocalizationPole inv(T(0))
  @test_throws ArgumentError GWF._cert_prime_type(9)
  @test_throws ArgumentError GWF._cert_prime_type(big(1)<<32)
  # A fractional answer must remain fractional. No integrality rounding is used.
  records = [Dict("accepted"=>true, "prime"=>101, "residues"=>[76])]
  recovery = GWF._cert_recover(records, big(4), [big(1)])
  @test recovery.values == [1//4]
  @test GWF._cert_recover(records, big(4), [big(1000)]) === nothing
  @test_throws ArgumentError GWF._cert_recover(vcat(records, records), big(4), [big(1)])
end

@testset "Certified arithmetic beyond the Peskine model" begin
  # Untwisted non-projective flag.
  X = Gr(2, 4)
  reference = small_quantum_multiplication(X; algorithm=:branches)
  certified = small_quantum_multiplication(X; arithmetic=:certified)
  @test certified.matrix == reference.matrix
  @test certified.metric == reference.metric

  # Ordinary split convex bundle.
  X = projective_space(5)
  E = line_bundle(X, 3)
  reference = small_quantum_multiplication(X; bundle=E, algorithm=:branches)
  certified = small_quantum_multiplication(X; bundle=E, arithmetic=:certified)
  @test certified.matrix == reference.matrix

  # An actual globally generated bundle with nonsplit graded pieces.
  # A general U3-dual section on Fl(1,3;4) gives P2, with degrees (d,0).
  X = flag_variety(4, [1, 3])
  L, B = tautological_bundles(X)
  twist = ConvexEulerTwist(dual(L) + dual(B);
    certificate="U3dual is globally generated; a general section has zero locus P2")
  degrees = [[0, 0], [1, 0]]
  reference = small_quantum_multiplication(
    X; bundle=twist, betas=degrees, algorithm=:branches
  )
  certified = small_quantum_multiplication(
    X; bundle=twist, betas=degrees, arithmetic=:certified
  )
  @test certified.matrix == reference.matrix
  @test size(certified.matrix) == (3, 3)
end

@testset "Exact flag integration and cyclic algebra" begin
  X = projective_space(1)
  value = equivariant_flag_integral(
    X; bundles=(line_bundle(X, 1), line_bundle(X, 0))
  ) do f, _
    @test all(iszero, f[2])
    [sum(f[1]), sum(f[1])+sum(f[2])]
  end
  @test value == [1, 1]
  A = [0 0; 1 0]
  algebra = cyclic_quantum_algebra(A, [1, 0]; metric=[0 1; 1 0])
  @test string(algebra.characteristic_polynomial)=="t^2"
  @test algebra.krylov_determinant==1
  @test length(algebra.blocks)==1
  @test only(algebra.blocks).dimension==2
  @test Matrix(algebra.multiplication_operators[2])==A
  @test_throws ArgumentError cyclic_quantum_algebra([1 0; 0 1], [1, 0])
end

@testset "Certified localization, resume, interruption and identity checks" begin
  X = projective_space(2)
  classes = [point_class(X)]
  mktempdir() do directory
    first_run = certified_pair_invariants(X, classes, [(1, 1)], [1]; checkpoint=directory)
    @test first_run.values == [1]
    @test first_run.certificate["new_prime_jobs"] > 0
    resume = certified_pair_invariants(X, classes, [(1, 1)], [1]; checkpoint=directory)
    @test resume.values == first_run.values
    @test resume.certificate["new_prime_jobs"] == 0
    @test resume.certificate["bound_reused"]
    @test_throws ArgumentError certified_pair_invariants(X, classes, [(1, 1)], [1];
      checkpoint=directory, torus_seeds=(0, 1))
    path = joinpath(directory, "beta-1", "bound.toml")
    record = TOML.parsefile(path)
    record["seconds"] += 1
    open(io -> TOML.print(io, record), path, "w")
    @test_throws ArgumentError certified_pair_invariants(
      X, classes, [(1, 1)], [1]; checkpoint=directory
    )
  end
  mktempdir() do directory
    interrupt =
      event -> event.phase == :certified_prime ? error("simulated interruption") : nothing
    @test_throws ErrorException certified_pair_invariants(X, classes, [(1, 1)], [1];
      checkpoint=directory, progress=interrupt)
    result = certified_pair_invariants(X, classes, [(1, 1)], [1]; checkpoint=directory)
    @test result.values==[1]
    @test result.certificate["bound_reused"]
    @test result.certificate["reused_prime_jobs"]>=1
  end
  @test_throws ArgumentError certified_pair_invariants(X, classes, [(1, 1)], [0])
  @test_throws ArgumentError certified_pair_invariants(
    X, classes, [(1, 1)], [1]; torus_seeds=(0, 0)
  )
  @test_throws ArgumentError certified_pair_invariants(
    X, classes, [(1, 1)], [1]; max_primes=0
  )
  @test_throws ArgumentError certified_pair_invariants(X, classes, [(1, 2)], [1])
  reference = small_quantum_multiplication(projective_space(1); algorithm=:branches)
  certified = small_quantum_multiplication(projective_space(1); arithmetic=:certified)
  @test certified.matrix == reference.matrix
  @test all(c -> c["status"]=="certified", certified.certification)
  @test_throws ArgumentError small_quantum_multiplication(X; arithmetic=:unknown)
end

if get(ENV, "GWFLAGS_PARALLEL_TESTS", "false")=="true"
  @testset "Certified worker isolation and parallel resume" begin
    before = Distributed.workers()
    mktempdir() do directory
      X = projective_space(2)
      result = certified_pair_invariants(X, [point_class(X)], [(1, 1)], [1];
        checkpoint=directory, workers=2)
      @test result.values==[1]
      @test Distributed.workers()==before
      resume = certified_pair_invariants(X, [point_class(X)], [(1, 1)], [1];
        checkpoint=directory, workers=1)
      @test resume.values==result.values
      @test resume.certificate["new_prime_jobs"]==0
      @test Distributed.workers()==before
    end
  end
end
