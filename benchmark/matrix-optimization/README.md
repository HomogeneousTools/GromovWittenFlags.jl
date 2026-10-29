# Local matrix optimization checks

Twenty additional fourfold models were checked on 29 October 2026 against
the unchanged engine at `e5753890ac1244eb2c9730d749df05e6b4d0c627`.
The candidate combined curve-budget pruning, smaller exact torus coordinates,
FLINT rational arithmetic in positive-degree quantum blocks, and reuse of
owned product temporaries. These kernels were tested first as process-local
overrides and then integrated into the package.

The sample covers twenty distinct geometric equivalence classes, ambient
dimensions 6–19, Picard ranks 2–5, and 24–630 fixed points. Eighteen bundles
have nonlinear summands and nine cases have disconnected zero loci. It is a
coverage sample, not an estimate of the remaining HPC workload.

## Exact references

`references.toml` retains every Novikov coefficient of the matrix, metric,
grading, basis keys, curve list, and the 101 ordinary and 101 regularized
period coefficients for each model. All twenty candidate results matched
exactly: 2,918 matrix entries, 1,020 nonzero Novikov coefficients, and 4,040
period coefficients per engine. A separate Python `Fraction` implementation
recovered the periods by graded triangular substitution and checked the
operator grading, metric, and self-adjointness.

Nineteen baseline matrices finished locally. The last local baseline hit a
one-hour process limit. Its exact reference is the existing completed HPC
output from the same pinned engine, verified locally against the prototype.
The archived record identifies this origin; its HPC timing is not included
in any local matrix speedup.

The package tests retain the original coloring iterator as a reference and
compare full decorated-orbit sets, multiplicities, and automorphism orders.
The six-vertex multigraded test covers 4,368 rooted/unrooted fixtures and
239,928 graph orbits. Independent AHU keys identify graphs and stabilizers.
Arithmetic tests also check that product reuse leaves cached factors intact.

## Timing method and results

Apple M4, Julia 1.13.0, PartialFlagVarieties 0.4.0, Nemo 0.54.2, one Julia
thread and one BLAS thread. Each engine/model pair used a separate process.
The first matrix call includes compilation. A second, warmed call was
omitted whenever the first call exceeded 300 seconds. Construction, period
extraction, comparison, and output lie outside the matrix timer.

| Matrix calls | Pairs | Median speedup | Summed-time speedup |
|---|---:|---:|---:|
| First call, including compilation | 19 | 2.54× | 5.48× |
| Warmed call | 15 | 3.52× | 3.81× |

The first-call range is 0.95–8.59×; all fifteen warmed comparisons improve
(1.11–8.87×). Two short first calls initially appeared slightly slower.
Two additional fresh processes per engine put their median speedups at
1.03× and 1.00×, consistent with parity when compilation dominates.

The heavy nineteenth case, `140600800.w2Bw142w0A`, changes from 3,399.93 to
395.85 seconds (8.59×). The final case, `12040A00.w1J4w0Zw0Zw06`, took
415.99 seconds for the prototype matrix and 424.98 seconds for its whole
worker. The unfinished baseline worker exceeded 3,600.14 seconds, giving
only a whole-process speedup bound of **greater than 8.47×**. It is excluded
from the nineteen finite matrix ratios.

Accepted timings exclude system suspension: direct macOS active-uptime clocks
are used where recorded; earlier elapsed timers are accepted only when process
clocks establish no suspension. An interrupted call without an active clock
was repeated. Julia allocation counters exclude some native allocations and
do not measure peak resident memory. These are local results; no optimized
HPC timings have been measured.

## Reproduction

Run one model per fresh process to measure compilation consistently:

```sh
julia --project=. --startup-file=no --threads=1 benchmark/matrix_optimization.jl 2230.7Z --warm
```

The script checks the complete output against the archived reference after
each timed matrix call. It writes ignored `local-results.toml` in this
directory. With no labels it checks all twenty cases in one process, which
shares compilation across models and therefore has a different timing protocol.
Use the same script and reference file in a separate baseline checkout to
repeat a paired comparison; the baseline may need substantially longer than
one hour for the last case. No HPC access is needed.
