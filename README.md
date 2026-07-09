# Bacteria WGS Pipeline

![Pipeline overview](pipeline.png)

One script for bacterial whole-genome sequencing analysis:
assembly → assembly QC → CheckM → Bakta annotation. `--mode` is required
and selects one of:

- **short** — Illumina paired-end only
  Trimmomatic → SPAdes → reformat.sh → QUAST → BBMap coverage
- **long_np** — Nanopore only
  Filtlong → Flye (`--nano-hq`/`--nano-raw`) → medaka → reformat.sh → QUAST
  → minimap2 (`map-ont`) + mosdepth coverage
- **long_pb** — PacBio HiFi only
  Filtlong → Flye (`--pacbio-hifi`) → reformat.sh → QUAST → minimap2
  (`map-hifi`) + mosdepth coverage
  (no polishing step — medaka is ONT-trained and doesn't apply to HiFi)
- **hybrid_np** — Illumina + Nanopore together
  Trimmomatic + Filtlong → SPAdes `--nanopore` → reformat.sh → QUAST →
  BBMap + minimap2 (`map-ont`) mapping, combined via pileup.sh
- **hybrid_pb** — Illumina + PacBio HiFi together
  Trimmomatic + Filtlong → SPAdes `-s` (HiFi fed in as an extra single-read
  library — SPAdes has no `--pacbio-hifi` flag; `--pacbio` is for noisy CLR
  reads) → reformat.sh → QUAST → BBMap + minimap2 (`map-hifi`) mapping,
  combined via pileup.sh

All modes converge on the same final assembly FASTA, so CheckM (batch, over
every sample) and Bakta annotation run identically regardless of mode.

## Requirements

The following tools must be on `PATH` (a working `environment.yml` that
installs all of them into one conda env is included):

- `trimmomatic`, `spades.py`, `bbmap.sh` — short / hybrid_*
- `filtlong`, `minimap2`, `samtools`, `flye`, `mosdepth` — long_* / hybrid_*
- `medaka_consensus` — long_np only
- `pileup.sh` (part of BBMap) — hybrid_* only
- `reformat.sh`, `quast.py`, `checkm`, `bakta` — every mode

```bash
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate bacteria_wgs
```

## Setup

`config.sh` holds machine-specific paths and is git-ignored, since it points
at things like your Trimmomatic adapter file and Bakta database — not
something to publish. Create your own from the template and edit it:

```bash
cp config.example.sh config.sh
```

Edit `config.sh` (or pass the equivalent flags instead) to set:

- `ADAPTER` — Trimmomatic adapter FASTA (short/hybrid_*)
- `BAKTA_DB` — Bakta database directory
- `TMP_DIR` — Bakta temp directory
- `THREADS`, `MIN_LENGTH` — defaults for every run
- `FILTLONG_MIN_LENGTH`, `FILTLONG_KEEP_PERCENT` — long-read filtering (long_*/hybrid_*)
- `FLYE_READ_TYPE` — `--nano-hq` / `--nano-raw` (long_np only)
- `MEDAKA_MODEL` — must match your basecaller/flow cell/kit; no safe default,
  the pipeline refuses to start in long_np mode until this is set
  (`medaka tools list_models` to see what's available)
- `FLYE_PACBIO_READ_TYPE` — `--pacbio-hifi` / `--pacbio-raw` (long_pb only)

## Running it from anywhere

To call the pipeline like a regular command instead of `cd`-ing into this
folder every time:

```bash
mkdir -p ~/tools/bacteria_wgs_pipeline
cp bacteria_wgs_pipeline.sh config.sh ~/tools/bacteria_wgs_pipeline/
chmod +x ~/tools/bacteria_wgs_pipeline/bacteria_wgs_pipeline.sh

conda activate bacteria_wgs
ln -s ~/tools/bacteria_wgs_pipeline/bacteria_wgs_pipeline.sh \
  "$CONDA_PREFIX/bin/bacteria_wgs_pipeline"
```

The script resolves symlinks back to its real directory, so it still finds
`config.sh` next to itself. From then on, whenever the `bacteria_wgs` env is
active, `bacteria_wgs_pipeline` works from any directory.

## Usage

`--mode` is required; the script refuses to run without it.

```bash
# Illumina only: every sample in raw_read/ (expects <sample>_1.fastq.gz / _2.fastq.gz)
./bacteria_wgs_pipeline.sh --mode short -r raw_read -o results

# Nanopore only: every sample in raw_long/ (expects <sample>.fastq.gz)
./bacteria_wgs_pipeline.sh --mode long_np -l raw_long -o results

# PacBio HiFi only
./bacteria_wgs_pipeline.sh --mode long_pb -l raw_long -o results

# Hybrid: matching samples in both directories
./bacteria_wgs_pipeline.sh --mode hybrid_np -r raw_read -l raw_long -o results
./bacteria_wgs_pipeline.sh --mode hybrid_pb -r raw_read -l raw_long -o results

# Process specific samples only
./bacteria_wgs_pipeline.sh --mode short -r raw_read -o results -s 18H3P11,18H3P12

# Re-run only the annotation step after fixing something upstream
./bacteria_wgs_pipeline.sh --mode short -r raw_read -o results --skip-checkm --force
```

Run `./bacteria_wgs_pipeline.sh -h` for the full option list.

## What it does

For every sample, depending on `--mode`:

- **short**: Trimmomatic (`01_QC/<sample>/`) → SPAdes `--isolate`
  (`02_assembly/<sample>/scaffolds.fasta`)
- **long_np**: Filtlong (`01_QC/<sample>/`) → Flye
  (`02_assembly/<sample>/assembly.fasta`) → medaka polishing
  (`02_assembly/<sample>/medaka/consensus.fasta`)
- **long_pb**: Filtlong → Flye `--pacbio-hifi`
  (`02_assembly/<sample>/assembly.fasta`, used directly — no polishing step)
- **hybrid_np**: Trimmomatic + Filtlong → SPAdes with `--nanopore`
  (`02_assembly/<sample>/scaffolds.fasta`)
- **hybrid_pb**: Trimmomatic + Filtlong → SPAdes with `-s` (HiFi as an
  extra single-read library) (`02_assembly/<sample>/scaffolds.fasta`)

Then, regardless of mode:

- **reformat.sh** drops contigs shorter than `MIN_LENGTH` (default 1000 bp)
  → `scaffolds_final.fasta`
- **QUAST** (`03_quast/<sample>/`) — assembly quality report, fed the
  relevant reads per mode (`--pe1/--pe2` for short reads, `--nanopore` for
  Nanopore, `--pacbio` for PacBio HiFi) so it also reports mapping-based
  metrics, not just bare assembly stats
- **Coverage** (`03_quast/<sample>/`), tool depends on mode:
  - short: `bbmap.sh covstats=` → `coverage_stats.txt`
  - long_np/long_pb: `minimap2` (`map-ont`/`map-hifi`) → sorted/indexed bam
    → `mosdepth` → `<sample>_coverage.mosdepth.summary.txt`
  - hybrid_*: short reads via `bbmap.sh` and long reads via `minimap2` are
    mapped to separate bams, then combined with `pileup.sh in=a.bam,b.bam`
    → `<sample>_total_hybrid_coverage.txt` (`samtools merge` doesn't play
    well with bbmap.sh vs minimap2 headers, so this sidesteps it)

Once every sample has a final assembly, they are all copied into
`04_checkm/bins/` and:

- **CheckM** runs once (`lineage_wf` + `qa`) over every sample's assembly
  together, matching how CheckM is normally used across a genome set.

Finally:

- **Bakta** annotates each sample's final assembly (`05_bakta/<sample>/`).

## Resuming / re-running

Each step is skipped if its expected output already exists, so a failed or
interrupted run can simply be re-invoked with the same arguments to continue
where it left off. Pass `--force` to re-run every step regardless. Logs for
every step are written to `00_logs/`.

## Notes

- CheckM is a batch step: it is most meaningful once you have multiple
  genomes to compare, so it only runs after all requested samples finish
  assembly. Use `--skip-checkm` to defer it (e.g. when adding a single new
  sample to a set already assessed).
