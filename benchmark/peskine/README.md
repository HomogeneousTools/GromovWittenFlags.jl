# Peskine reproducibility milestone

The portable example in `examples/peskine` is a demanding regression case for
the certified backend. Keep case data and scientific assertions outside `src/`.
Its coefficient fixture is consulted by the post-computation verifier, never
by the numerical engine.

On 17 November 2026 the complete prototype was regenerated on Aion with Julia
1.11.4, Nemo 0.54.2, PartialFlagVarieties 0.4.0 and Semisimple 1.3.0.
It computed all five bounds and 158 new prime jobs at two torus points, with
no reused computational inputs. Every coefficient matched the original exact
certificate. The trace was -633; the inherited and vanishing scalar actions
were -14 and -15. The resulting primary factors and the stated geometric
inputs yield coarse Hodge-atom dimensions 24,35,1.

The four-worker computation took 1954.963 seconds including package loading,
worker startup, fresh classical initialization, all certified degrees, native
trace and algebra, and output. Dependency installation was excluded. Including
the final comparison and polynomial export, the Slurm job took 32m50s.
Peak sampled sum of process RSS was 21.433 GiB at 0.1-second intervals; shared
pages can be counted more than once and brief peaks can be missed. Slurm's
separate MaxRSS measurement was 24.071 GiB. The memory request caused the
scheduler to allocate more CPUs than the four numerical workers used.

A new coordinator with two workers reproduced the answer from copied
checkpoints in 126.220 seconds, reusing all bounds and generating no new prime
jobs. That is a replay timing, not a second fresh-computation timing.

`milestone.toml` records public numerical facts and prototype engine hashes.
The source was subsequently formatted and documented for release. The timing
is a historical prototype measurement, not a new benchmark of the exact
release tree. Its source-bound checkpoint identities are retained with the
original separate scientific archive, not mixed with release checkpoints.
No account-specific Slurm launchers, private paths or runtime depots are part
of this package example.
