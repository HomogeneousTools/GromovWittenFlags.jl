# Contributing

Use the conventions of PartialFlagVarieties.jl: two-space indentation,
JuliaFormatter's Blue style, focused source files, docstrings on public and
private functions and types, and exact mathematical regression tests.

## Structure

`src/GromovWittenFlags.jl` declares explicit dependency imports and exports.
`Schubert.jl` handles flag root data and curve coordinates; `Bundles.jl` handles
Levi weights and twisting data. Tree enumeration is independent of localization.
`Localization.jl` owns numeric restrictions and graph factors; `QuantumBlocks.jl`
shares work between pairs of insertions. `Invariants.jl` validates public inputs
and chooses the scalar or pair-block calculation. `Novikov.jl` and
`ExactLinearAlgebra.jl` supply the exact algebra used by `Quantum.jl`.
`ExactSums.jl` provides balanced rational accumulation for graph contributions.
Only public APIs of PartialFlagVarieties and Semisimple are imported.

Use ordinary Julia dispatch and concrete arithmetic buffers. Put recursive
search state in explicit worker arguments; avoid recursive local functions or
captured variables that are reassigned elsewhere in the enclosing scope.
Use `!` for helpers that overwrite caller-owned workspaces, and document the
buffers and return value. A memoized reader must identify its cache effects.
Keep deliberate heterogeneous storage at a documented boundary, such as the
variety cache or TOML reporting, and explain why it is needed there.

## Tests and documentation

Run `julia --project=. -e 'using Pkg; Pkg.test()'` before submitting a change.
Tests should include independently known mathematical answers, input guards,
and symmetry multiplicities. Avoid comparisons solely against the Python port.
The full degree-five Küchle matrices are enabled with
`GWFLAGS_EXTENDED_TESTS=true`; manual CI dispatch runs them separately. The smaller full
fourfold matrices and Küchle checks through degree three belong to the default suite.

Build the strict Documenter manual using the commands in the README.
Documenter checks every docstring, including the private reference.
The suite also checks that private functions and types have documentation.
Explain mathematical signs, symmetry factors, borrowed buffers, and pruning
conditions inline where they are used. Examples must execute.
Docstrings should identify return types and index conventions where ambiguous.
Public examples belong in executable doctest blocks. Private contracts should
state preconditions and storage lifetimes rather than repeat the function name.
For a mutating helper, distinguish its returned value from the buffers it fills.
Mathematical regressions should compare exact objects and coefficients; keep
text rendering assertions in focused display tests.
Changes to curve or basis conventions require an explicit migration note.

Run `julia --project=. benchmark/benchmarks.jl` to compare warmed timings and
allocations. See [benchmark/README.md](benchmark/README.md) for the workloads
and measurement method. Preserve exact output checks when optimizing a sum;
two implementations of the same recurrence are weaker evidence than a small
brute-force oracle or an independently known invariant.

## Formatting

`julia devtools/format.jl` formats source, tests, examples, and developer tools.
`julia devtools/format.jl --check` checks without changing files.
The formatter version is pinned in that script; local checks and CI use it.

Benchmark scripts define their tools in dedicated modules and run only through
their CLI guards. The GKMtools comparison has its own environment and tests;
see [its instructions](benchmark/gkmtools/README.md). External comparison
dependencies must not enter the production package or its ordinary test suite.

## Mathematical review

Parts of this implementation were developed with language-model assistance.
Human review is essential before merging mathematical changes. Preserve source
provenance and third-party notices. A dimension check and agreement at two
torus points provide useful checks, but do not replace review of the formulas.
