# Comparing completed HPC cases locally

`hpc_comparison.jl` measures a complete quantum multiplication matrix with either
the current package or an unchanged archived `GromovWitten.jl` engine. It uses
the same installed geometry dependencies for both engines. The old file is
included in a separate module with the imports formerly supplied by its parent
package; its bytes must match the SHA-256 recorded in the input.

Prepare one local TOML input per completed fourfold case:

```toml
label = "50.F" # Quartic fourfold in projective five-space.
ambient_factors = [[1, 6]]
anticanonical_coordinates = [2]
betas = [[0], [1], [2]]
old_source = "/absolute/path/to/the/archived/GromovWitten.jl"
old_source_sha256 = "the SHA-256 recorded by the HPC run"
new_source_commit = "the tested package commit"
```

The geometry fields above give a small smoke test. Fill in the archive path,
hash and tested commit. For a completed HPC workload, copy the actual factors,
coordinates and complete beta list from the case description and saved output.
The worker independently enumerates the curve classes and rejects a different
list. No degree cap is applied.

Run each backend in a fresh process, serially, on the same otherwise idle host:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p julia --project=. --startup-file=no \
  --threads=1 --gcthreads=1 --heap-size-hint=1G \
  benchmark/hpc_comparison.jl old INPUT.toml OLD.toml
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p julia --project=. --startup-file=no \
  --threads=1 --gcthreads=1 --heap-size-hint=1G \
  benchmark/hpc_comparison.jl new INPUT.toml NEW.toml
```

The dependency manifest must be identical. Both geometry cache budgets are
128 MiB and BLAS uses one thread. Alternate backend order between cases. The
`time -p` option reports whole-process time. On macOS, `time -l` also reports
peak resident memory when the environment permits access to `kern.clockrate`;
restricted sandboxes can reject that query even after the worker succeeds.

Each process constructs the target, warms the classical term and the first
nonzero curve, performs a full garbage collection, then measures **one complete
matrix computation**. Construction, warm-up, and the period recurrence are
reported separately. Remaining compilation is included in matrix time and
recorded explicitly; these are single samples, not medians or fully warmed
measurements. Allocated bytes measure cumulative allocation, not resident memory.

Set `warmup = false` in the input TOML to measure the first full matrix call
in a fresh process. Geometry construction and automatic curve enumeration still
precede the matrix timer. Record these calls separately from measurements after
the short warm-up: another implementation can populate substantial numerical
caches during that warm-up. First-call timings include matrix compilation.

The output includes every exact Novikov coefficient, the metric, basis keys,
codimensions, grading, curve classes, and both ordinary and regularized period
coefficients through order 100. Compare all of these against the saved HPC
output and between engines. Polynomial entries are strings of sorted
`exponent=coefficient` pairs separated by semicolons; the empty string is zero.
Rationals use `numerator/denominator`, with denominator one omitted. The HPC
JSON grading may need its `//` spelling normalized before comparison.

The worker checks the connection grading before applying the same fourfold
period recurrence as the HPC runner. This is a consistency check on the ambient
quantum connection; it does not reconstruct primitive cohomology. Completed
results and the preceding warm-up checkpoint are written atomically.

Compare local old/new times to measure engine changes. The original HPC times
identify the workloads; differences between HPC and local times also reflect
hardware, Julia versions, contention, and compilation. Three cases cannot
establish a campaign-wide speedup.
