# bacteria-wgs-pipeline

Bacterial WGS pipeline: bash script (`bacteria_wgs_pipeline.sh`) plus a
Nextflow DSL2 port (`nextflow/`). Five modes: short, long_np, long_pb,
hybrid_np, hybrid_pb. See `README.md` for details.

## Code philosophy (ponytail-inspired)

Before adding code, work down this ladder and stop at the first "yes":

1. Does this need to exist at all? (YAGNI)
2. Does it already exist in this codebase?
3. Does the standard library / already-installed tool do this?
4. Can it be one line?
5. Only then: the minimum that actually works.

Prefer deletion over addition, boring over clever, fewest files possible.
When two options are the same size, pick the more robust one over the
merely simpler-looking one.

Non-negotiable regardless of the above: understand the problem before
coding, validate input at trust boundaries, handle errors that could lose
data, and verify non-trivial logic actually runs (a manual test, a stub
run, or equivalent) before calling it done.
