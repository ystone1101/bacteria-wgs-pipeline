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

## Getting started

### 1. Clone the repository

```bash
git clone https://github.com/ystone1101/bacteria-wgs-pipeline.git
cd bacteria-wgs-pipeline
```

### 2. Create the conda environment

Requires [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/)
(mamba solves much faster — recommended). `environment.yml` installs every
tool the pipeline needs into one env called `bacteria_wgs`:

```bash
mamba env create -f environment.yml   # or: conda env create -f environment.yml
```

### 3. Activate it and sanity-check the tools

```bash
conda activate bacteria_wgs

trimmomatic -version
spades.py --version
flye --version
medaka --version
quast.py --version
checkm -h > /dev/null && echo "checkm OK"
bakta --version
```

You'll need to reactivate this env (`conda activate bacteria_wgs`) in every
new terminal session before running the pipeline.

### 4. Set up `config.sh`

`config.sh` holds machine-specific paths (adapter file, Bakta database,
etc.) and is git-ignored — it's meant to stay local, not get committed.
Create your own from the template:

```bash
cp config.example.sh config.sh
```

Then edit `config.sh` and fill in:

| Variable | Used by | Notes |
|---|---|---|
| `THREADS` | every mode | default thread count for every tool |
| `ADAPTER` | short / hybrid_* | Trimmomatic adapter FASTA, absolute path |
| `BAKTA_DB` | every mode | Bakta database directory, absolute path |
| `TMP_DIR` | every mode | Bakta temp directory, absolute path |
| `MIN_LENGTH` | every mode | min contig length kept after assembly |
| `FILTLONG_MIN_LENGTH`, `FILTLONG_KEEP_PERCENT` | long_*/hybrid_* | long-read filtering |
| `FLYE_READ_TYPE` | long_np | `--nano-hq` or `--nano-raw` |
| `MEDAKA_MODEL` | long_np | **no safe default** — must match your basecaller/flow cell/kit, or the pipeline refuses to start; run `medaka tools list_models` to see what's available |
| `FLYE_PACBIO_READ_TYPE` | long_pb | `--pacbio-hifi` or `--pacbio-raw` |

Any of these can also be overridden per-run with a command-line flag
instead of editing `config.sh` (see `-h`).

### 5. (Optional) Put the pipeline on `PATH`

So you can call it as a plain command from any directory, instead of
`cd`-ing into this folder every time:

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
active, `bacteria_wgs_pipeline` works from any directory. (The examples
below use `./bacteria_wgs_pipeline.sh`; swap in `bacteria_wgs_pipeline` if
you did this step.)

### 6. Lay out your reads

Sample names come from the read filenames — pick a folder per read type
and name files by sample:

- Short reads (Illumina): `<sample>_1.fastq.gz`, `<sample>_2.fastq.gz`
- Long reads (Nanopore/PacBio): `<sample>.fastq.gz` (one file per sample)

For `hybrid_np`/`hybrid_pb`, the short-read and long-read directories must
use the **same sample name** for the same sample — that's how the script
matches an Illumina pair to its long reads:

```
raw_read/16N2L7_1.fastq.gz   raw_read/16N2L7_2.fastq.gz
raw_long/16N2L7.fastq.gz
```

### 7. Run it

`--mode` is required; the script refuses to run without it.

```bash
# Illumina only: every sample in raw_read/
./bacteria_wgs_pipeline.sh --mode short -r raw_read -o results

# Nanopore only: every sample in raw_long/
./bacteria_wgs_pipeline.sh --mode long_np -l raw_long -o results

# PacBio HiFi only
./bacteria_wgs_pipeline.sh --mode long_pb -l raw_long -o results

# Hybrid: matching samples in both directories
./bacteria_wgs_pipeline.sh --mode hybrid_np -r raw_read -l raw_long -o results
./bacteria_wgs_pipeline.sh --mode hybrid_pb -r raw_read -l raw_long -o results

# Process specific samples only
./bacteria_wgs_pipeline.sh --mode short -r raw_read -o results -s 18H3P11,18H3P12
```

Run `./bacteria_wgs_pipeline.sh -h` for the full option list.

### 8. Check the results

Everything lands under the `-o` directory, numbered so it sorts in
pipeline order:

```
results/
├── 00_logs/         one log file per sample per step
├── 01_QC/<sample>/  trimmed/filtered reads
├── 02_assembly/<sample>/
├── 03_quast/<sample>/  QUAST report + coverage stats
├── 04_checkm/       bins/ (all samples' assemblies) + output/ (batch results)
└── 05_bakta/<sample>/  annotation (gff3, faa, gbff, ...)
```

Progress prints with colored `START`/`DONE`/`FAIL`/`SKIP` status and a
summary at the end; colors auto-disable when output isn't a terminal (e.g.
redirected to a log file).

### 9. If a run fails partway through

Just re-run the same command. Each step is skipped if its expected output
already exists, so you resume right after the failure instead of starting
over. Pass `--force` to re-run every step regardless.

## What each mode does, in detail

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

## Notes

- CheckM is a batch step: it is most meaningful once you have multiple
  genomes to compare, so it only runs after all requested samples finish
  assembly. Use `--skip-checkm` to defer it (e.g. when adding a single new
  sample to a set already assessed).
- `bakta<=1.11.4` crashes late in the annotation step (`AttributeError:
  'str' object has no attribute 'decode'`) if paired with `pyhmmer>=0.12`
  — `environment.yml` pins `pyhmmer<0.12` to avoid this.
