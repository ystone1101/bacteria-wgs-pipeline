#!/usr/bin/env bash
# Default settings for bacteria_wgs_pipeline.sh
# Edit the paths below to match your environment, or override any of them
# on the command line (see ./bacteria_wgs_pipeline.sh -h).

# Number of threads to use for every tool
THREADS=6

# Trimmomatic adapter FASTA (ILLUMINACLIP)
ADAPTER="/path/to/trimmomatic_adapters/TruSeq3-PE-2-GGGGG.fa"

# Bakta database directory
BAKTA_DB="/path/to/bakta_db/db"

# Temp directory used by Bakta
TMP_DIR="/tmp/bakta_tmp"

# Minimum contig length kept after assembly (reformat.sh minlength=)
MIN_LENGTH=1000

# --- long_np/long_pb/hybrid_* mode only ---

# Filtlong read-level filtering (long-read QC)
FILTLONG_MIN_LENGTH=1000
FILTLONG_KEEP_PERCENT=90

# --- long_np (Nanopore) only ---

# Flye read type: --nano-hq for recent high-accuracy basecalls,
# --nano-raw for older/fast basecalling models
FLYE_READ_TYPE="--nano-hq"

# If a run fails with Flye's "No disjointigs were assembled" error (seen
# on very high per-sample coverage, e.g. deep PromethION runs), pass
# -g <genome size> --asm-coverage <N> on the command line (see -h) --
# these aren't config.sh settings since genome size varies per sample.

# medaka polishing model - must match your basecaller/flow cell/kit,
# e.g. "r1041_e82_400bps_sup_v5.0.0". No safe default; run
# `medaka tools list_models` to see what's available.
MEDAKA_MODEL=""

# --- long_pb (PacBio) only ---

# Flye read type: --pacbio-hifi for HiFi reads, --pacbio-raw for older CLR reads
FLYE_PACBIO_READ_TYPE="--pacbio-hifi"
