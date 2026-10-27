# GKMtools comparison

This isolated environment compares exact first-Chern quantum multiplication
on Fano fourfolds. GKMtools and Oscar are benchmark dependencies only.
The production package continues to depend on PartialFlagVarieties and
Semisimple. See [RESULTS.md](RESULTS.md) for the measured results and decision.

## Reproduction

From the package root, using Julia 1.13:

`Manifest-v1.13.toml` records the complete dependency resolution used for the
measurements. Other Julia versions resolve from `Project.toml` instead.

```sh
julia --project=benchmark/gkmtools -e 'using Pkg; Pkg.instantiate()'
JULIA_NUM_THREADS=1 julia --compiled-modules=existing --startup-file=no \
  --project=benchmark/gkmtools benchmark/gkmtools/fourfolds.jl \
  cubic Gushel-Mukai quartic quadric-cubic three-quadrics
```

The recorded run used `--compiled-modules=existing` because normal Oscar
precompilation failed in the installed Julia 1.13 environment, including
Polymake and its Hecke extension. This option uses existing images and loads
other modules from source. Package loading is outside all measurements.
The results do not measure startup or compare Julia versions.

To compare the prepared integration adapter, retaining its classical pairing:

```sh
GWFLAGS_BACKENDS=GromovWittenFlags,GKMtools-prepared JULIA_NUM_THREADS=1 \
  julia --compiled-modules=existing --startup-file=no \
  --project=benchmark/gkmtools benchmark/gkmtools/fourfolds.jl \
  cubic Gushel-Mukai quartic quadric-cubic three-quadrics c5 b8
```

For a product with higher Picard rank, inspect each curve class separately:

```sh
GWFLAGS_TOTAL_DEGREE=2 GWFLAGS_SPLIT_DEGREES=true JULIA_NUM_THREADS=1 \
  julia --compiled-modules=existing --startup-file=no \
  --project=benchmark/gkmtools benchmark/gkmtools/fourfolds.jl K3-29
```

`GWFLAGS_BENCH_SAMPLES` sets a positive odd sample count, default three.
`GWFLAGS_BENCH_OUTPUT` selects the TOML output, default `results-local.toml`
next to the script. Each completed backend is saved immediately; an interrupted
backend is recorded without inventing a median. A whole-process termination
cannot save the interrupted row.

`GWFLAGS_KUECHLE_DEGREE` caps c5 and b8, default three. Their complete operators
have degree-five terms. `GWFLAGS_TOTAL_DEGREE` optionally caps every example
by the sum of its ambient curve coordinates. Without it, the six smaller
examples use their complete automatic degree lists. Setting
`GWFLAGS_GKM_BASIS=schur` selects the less sparse Schur-polynomial lift control.
`GWFLAGS_GKM_PREPARED=true` adds the prepared adapter to both ordinary backends;
`GWFLAGS_BACKENDS` can select any subset explicitly.

Run the independent root-data compatibility checks with:

```sh
julia --compiled-modules=existing --startup-file=no \
  --project=benchmark/gkmtools benchmark/gkmtools/weyl_checks.jl
```

## What is compared

`Models.jl` constructs each geometry independently with the native constructors
of each library. GKMtools uses its fast type-A permutation graphs and its native
tautological, quotient, exterior, symmetric, and line bundles. It does not use
PartialFlagVarieties to obtain its graph or bundle fibers.

We transport the selected ambient quotient basis by reduced words and the
Grassmannian partition convention. The default GKMtools restrictions are its
own Billey recurrence. Its public Schubert constructor requires Weyl-backed
vertices, whereas its fast type-A constructor stores permutation tuples.
Two private input-preparation helpers from the pinned revision bridge these
representations. They are not used in timed integration. Projective product
restrictions use the explicit product of their factor Schubert restrictions.
The alternative Schur lifts represent the same ordinary classes.

The independently computed twisted classical intersection matrices must agree
exactly before timing. Primitive curves are matched through a simple-root edge
or the two factor hyperplane degrees. Every individual Novikov coefficient
matrix, including degree zero, must agree after warm-up and after **every**
timed sample. Distinct curve classes are never combined by evaluating at q=1.

Three computational paths are distinguished:

- `GromovWittenFlags`: the public `small_quantum_multiplication` API, including
  quotient-basis selection, classical pairing, and both exact torus evaluations.
- `GKMtools`: its public `twisted_matrix_small_quantum_product`, called once
  for each requested curve class with the common basis. Its nonequivariant
  mode specializes to exact FLINT rational numbers. It rebuilds the symbolic
  classical pairing on each call.
- `GKMtools-prepared`: the small adapter in `PreparedOperator.jl`. It retains
  a native GKMtools classical matrix and inverse pairing, selects the
  dimension-compatible pairs, and calls its public `gromov_witten_nomarks`
  with `virtual_zero_section`. This gives GKMtools a deliberately favorable
  comparison, excluding its separately recorded classical preparation.
  GromovWittenFlags still recomputes its classical data inside the timed API.

All paths use the same requested degrees and common basis. Output conversion
and equality checks are outside the timed sections. There is one warm-up,
then three serial samples with a full collection before each. The file records
all elapsed times, Julia-reported allocated bytes, medians, package versions,
and the pinned GKMtools revision. Allocated bytes are not peak resident memory
and may omit allocations inside native libraries.

GKMtools uses a new rational torus specialization per positive curve class.
The harness fixes the random seed to 20261027 for warm-up and increments it
for each sample. GromovWittenFlags' built-in two-point agreement check remains
enabled. The extra cross-library comparisons are outside timing; they are not
a claim that GKMtools itself checks two torus points per result.

Combined `setup_seconds` includes both native models, basis transport, metric
verification, and any compilation encountered there. It cannot be interpreted
as either library's independent startup or setup cost. For the prepared path,
`gkm_preparation_seconds` records its additional native classical preparation.

## Why keep this separate

GKMtools has a broader remit, including generic GKM spaces, equivariant output,
descendants, positive-genus routines, and an untwisted quantum Schubert engine.
These tests address the convex Euler-twisted genus-zero fourfold workload here.
They do not rank those other features or justify replacing the Semisimple
root-data layer in PartialFlagVarieties.
