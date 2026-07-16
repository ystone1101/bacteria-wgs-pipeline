# bacteria-wgs-pipeline (Nextflow)

A personal Nextflow DSL2 port of [`bacteria_wgs_pipeline.sh`](../README.md). Same
five modes, same tools, same command lines under the hood -- just run through
Nextflow instead of a bash loop, so steps for different samples can run in
parallel and a crashed run resumes with `-resume` instead of marker files.

> **Execution status:** `short` and `long_np` have both run end-to-end on
> real data through every step -- `short`: Trimmomatic → SPAdes →
> reformat.sh → QUAST → BBMap coverage → CheckM → Bakta; `long_np`:
> Filtlong → Flye → medaka → reformat.sh → QUAST → minimap2/mosdepth
> coverage → CheckM → Bakta -- all completing successfully. `long_pb`,
> `hybrid_np`, `hybrid_pb` are not yet execution-tested through Nextflow
> (their bash equivalents all work). Do a small dry run on one sample before
> trusting an untested mode on a real dataset, and compare its output
> against the bash version's if anything looks off.
>
> Getting `short` mode running took a few real fixes along the way, all
> already applied here: Nextflow's newer versions (26.x) use a stricter
> parser that rejects `switch/case` and top-level statements outside a
> process/workflow (use `if/else` and keep everything inside `workflow {}`);
> and if your Nextflow engine needs an older JDK than your tools' conda env
> (`NXF_JAVA_HOME` set to boot Nextflow itself), that JDK's `bin/` leaking
> onto every task's `PATH` can break Java-based tools like Trimmomatic --
> `nextflow.config`'s `process.beforeScript` here clears that before each
> task runs. If you hit `Unable to parse config file ... Unsupported class
> file major version`, pin an older Nextflow via `export NXF_VER=24.10.5`.

## Modes

Same as the bash pipeline -- `short`, `long_np`, `long_pb`, `hybrid_np`,
`hybrid_pb`. See the main [README](../README.md) for what each mode does and
why (SPAdes flags, medaka/PacBio, etc.).

## Requirements

- Nextflow >= 23.04
- The same tools as the bash pipeline, available either via the shared conda
  env (`-profile conda`, reuses `../environment.yml`) or Docker/Singularity
  (`-profile docker` / `-profile singularity` -- no container images are
  pinned yet, you'll need to supply your own via `process.container` if you
  go this route)
- A Bakta database and (unless `--skip_checkm`) a CheckM reference database
  with `checkm data setRoot` already run once for your user -- see the main
  README's "Set up reference databases" step; this is identical for both the
  bash and Nextflow versions.

## Samplesheet

Instead of `-r`/`-l` directories, the Nextflow version takes one CSV via
`--input`:

```csv
sample,fastq_1,fastq_2,long_reads
S1,/full/path/to/raw_read/S1_1.fastq.gz,/full/path/to/raw_read/S1_2.fastq.gz,/full/path/to/raw_long/S1.fastq.gz
S2,/full/path/to/raw_read/S2_1.fastq.gz,/full/path/to/raw_read/S2_2.fastq.gz,/full/path/to/raw_long/S2.fastq.gz
```

See `assets/samplesheet_example.csv`. Rules:

- Paths must be readable from wherever Nextflow runs (absolute paths are
  safest).
- `short` mode: leave `long_reads` blank.
- `long_np` / `long_pb` mode: leave `fastq_1`/`fastq_2` blank.
- `hybrid_np` / `hybrid_pb` mode: fill in all three columns.

## Running it

```bash
# short (Illumina only)
nextflow run main.nf -profile conda \
  --mode short --input samplesheet.csv \
  --adapter /abs/path/adapters.fa --bakta_db /abs/path/bakta_db

# long_np (Nanopore only)
nextflow run main.nf -profile conda \
  --mode long_np --input samplesheet.csv \
  --medaka_model r941_prom_hac_g507 --bakta_db /abs/path/bakta_db

# long_pb (PacBio HiFi only, no medaka)
nextflow run main.nf -profile conda \
  --mode long_pb --input samplesheet.csv --bakta_db /abs/path/bakta_db

# hybrid_np / hybrid_pb (Illumina + long reads)
nextflow run main.nf -profile conda \
  --mode hybrid_np --input samplesheet.csv \
  --adapter /abs/path/adapters.fa --bakta_db /abs/path/bakta_db
```

Useful flags (all set in `nextflow.config`'s `params` block, override with
`--flag value`): `--outdir`, `--threads`, `--min_length`,
`--filtlong_min_length`, `--filtlong_keep_percent`, `--flye_read_type`,
`--flye_pacbio_read_type`, `--skip_checkm`, `--skip_bakta`.

To resume after a crash or an added sample, add `-resume` -- Nextflow's own
task-hash caching replaces the bash script's marker-file skip logic.

## Output layout

Same numbering as the bash pipeline, under `--outdir` (default `results`):
`00_logs/`, `01_QC/`, `02_assembly/`, `03_quast/`, `04_checkm/`,
`05_bakta/`, one subdirectory per sample (CheckM's output is batched
across all samples, same as the bash version).

`00_logs/` holds each tool's full stdout/stderr (named like the bash
version's, e.g. `<sample>_coverage.log`, `<sample>_medaka.log`) for every
step that doesn't already write its own log file into its output
directory -- SPAdes and Flye do that natively (`spades.log`/`flye.log`
alongside the assembly under `02_assembly/<sample>/`), so those aren't
duplicated into `00_logs/`.

## Differences from the bash version

- Samplesheet CSV instead of `-r`/`-l` directories with filename-convention
  sample discovery.
- Parallel execution and `-resume` instead of sequential steps with
  marker-file skipping.
- No `-s`/`-t`/`-a`/`-d`/`-m`/`--tmp-dir`/`--force` short flags -- use the
  `--param value` form throughout.
- Bakta's `--tmp-dir` is per-task (inside the Nextflow work dir) rather than
  one shared directory, so parallel samples can't collide.
