# Performance measurements

Run from the package checkout:

```sh
julia --project=. benchmark/benchmarks.jl
GWFLAGS_BENCH_SAMPLES=7 julia --project=. benchmark/benchmarks.jl
```

The script warms each workload once, then reports the median elapsed time
and median allocated bytes from an odd number of samples (three by default).
It runs a full garbage collection before each sample. Compilation is excluded;
the first process invocation can still spend time precompiling dependencies.
Use the same Julia version, dependency manifest, and otherwise idle machine
for comparisons. These measurements describe the selected workloads; they
are not promises about larger curve degrees or other targets.

## Review comparison

Measured on the same macOS host with Julia 1.13.0, sequentially, using three
warmed samples per workload. Both checkouts used the same installed dependency
versions. The baseline is the package immediately after extraction, before
the subsequent review. Times are seconds; allocations are bytes.
The dependency versions were PartialFlagVarieties 0.4.0, Semisimple 1.3.0,
StaticArrays 1.9.22, and Combinatorics 1.1.0.

| Workload | Before | After | Speedup | Before allocations | After allocations |
|:--|--:|--:|--:|--:|--:|
| Five-point degree-two invariant on the projective plane | 0.108921 | 0.000269 | 405.0× | 50,191,504 | 450,800 |
| Küchle c5 classical matrix | 0.313578 | 0.008272 | 37.9× | 457,711,224 | 18,537,976 |
| Küchle c5 matrix through curve degree two | 0.622733 | 0.219086 | 2.8× | 793,838,656 | 289,665,792 |
| Küchle c5 matrix through curve degree three | 4.546747 | 2.524987 | 1.8× | 4,558,133,848 | 2,004,059,240 |
| Quartic fourfold, full automatic degree list | 0.012153 | 0.004031 | 3.0× | 18,178,064 | 8,061,976 |
| K3-29 fourfold, full automatic degree list | 0.747074 | 0.202365 | 3.7× | 1,286,165,736 | 310,862,664 |

The conic workload benefits most from factoring labeled mark placements into
vertex sums. Shared classical restrictions, incremental row selection, and
reused LU factorizations account for the classical improvement. Covered-edge
and vertex caches, normalized Euler factors, and shared pair sums reduce the
positive-degree work. Balanced streaming sums limit intermediate denominators.
All calculations still use exact `Rational{BigInt}`
arithmetic and compare the full result at two torus points.

The test suite verifies the known outputs for these geometries. Independent
small subword and marked-orbit enumerations check the optimized kernels.
The script is intentionally dependency-free beyond the package dependencies;
for detailed profiling, wrap individual workloads in Julia's profiler.
