#!/usr/bin/env bash
#
# bacteria_wgs_pipeline.sh
#
# Bacterial WGS pipeline. One script, --mode is required:
#
#   short      Illumina paired-end only
#              Trimmomatic -> SPAdes -> reformat.sh -> QUAST -> BBMap
#
#   long_np    Nanopore long reads only
#              Filtlong -> Flye -> medaka -> reformat.sh -> QUAST
#   long_pb    PacBio HiFi long reads only
#              Filtlong -> Flye (--pacbio-hifi) -> reformat.sh -> QUAST
#              (no medaka: medaka is an ONT-trained polisher, not used for
#              already-accurate HiFi reads)
#
#   hybrid_np  Illumina + Nanopore together
#   hybrid_pb  Illumina + PacBio HiFi together
#              Trimmomatic + Filtlong -> SPAdes (--nanopore, or -s for HiFi)
#              -> reformat.sh -> QUAST -> BBMap
#
# All modes converge on the same final assembly FASTA, so CheckM (batch,
# over every sample) and Bakta annotation run identically regardless of mode.
#
# Usage:
#   ./bacteria_wgs_pipeline.sh --mode short     -r raw_read -o results
#   ./bacteria_wgs_pipeline.sh --mode long_np   -l raw_long -o results
#   ./bacteria_wgs_pipeline.sh --mode hybrid_pb -r raw_read -l raw_long -o results
#
# Run ./bacteria_wgs_pipeline.sh -h for the full option list.

set -euo pipefail

# Resolve the script's real directory, following symlinks, so that
# config.sh is still found when this script is invoked via a symlink
# on PATH (e.g. ~/bin/bacteria_wgs_pipeline -> .../bacteria_wgs_pipeline.sh).
resolve_script_dir() {
  local source="${BASH_SOURCE[0]}"
  while [[ -h "$source" ]]; do
    local dir
    dir=$(cd -P "$(dirname "$source")" && pwd)
    source=$(readlink "$source")
    [[ "$source" != /* ]] && source="$dir/$source"
  done
  cd -P "$(dirname "$source")" && pwd
}
SCRIPT_DIR=$(resolve_script_dir)

# ---------------------------------------------------------------------------
# Defaults (overridable via config.sh and/or command-line flags)
# ---------------------------------------------------------------------------
MODE=""
RAW_DIR="raw_read"
LONG_DIR="raw_long"
OUT_DIR="results"
SAMPLES=""
SKIP_CHECKM=0
SKIP_BAKTA=0
FORCE=0
GENOME_SIZE=""
ASM_COVERAGE=""

if [[ -f "$SCRIPT_DIR/config.sh" ]]; then
  # shellcheck source=config.sh
  source "$SCRIPT_DIR/config.sh"
elif [[ -f "$SCRIPT_DIR/config.example.sh" ]]; then
  echo "WARNING: config.sh not found next to the script; using placeholder" \
       "paths from config.example.sh. Run: cp config.example.sh config.sh" \
       "and edit it." >&2
  source "$SCRIPT_DIR/config.example.sh"
else
  echo "ERROR: neither config.sh nor config.example.sh found next to the script." >&2
  exit 1
fi

usage() {
  cat <<EOF
Usage: $(basename "$0") --mode {short,long_np,long_pb,hybrid_np,hybrid_pb} [options]

  --mode MODE   REQUIRED. One of:
                  short      Illumina only
                  long_np    Nanopore only
                  long_pb    PacBio HiFi only
                  hybrid_np  Illumina + Nanopore
                  hybrid_pb  Illumina + PacBio HiFi
  -r DIR        Short-read directory: <sample>_1.fastq.gz/_2.fastq.gz,
                <sample>_R1.fastq.gz/_R2.fastq.gz, or Illumina's own
                <sample>_S#_L#_R1_001.fastq.gz/_R2_001.fastq.gz (single
                lane only -- concatenate multi-lane runs first)
                (short/hybrid_* modes; default: ${RAW_DIR})
  -l DIR        Long-read directory, one <sample>.fastq.gz per sample
                (long_*/hybrid_* modes; default: ${LONG_DIR})
  -o DIR        Output directory for all pipeline results (default: ${OUT_DIR})
  -s LIST       Comma-separated sample IDs to process
                (default: auto-detect from the relevant read directory)
  -t THREADS    Threads for every tool (default: ${THREADS})
  -a FILE       Trimmomatic adapter FASTA (default: ${ADAPTER})
  -d DIR        Bakta database directory (default: ${BAKTA_DB})
  -m INT        Minimum contig length kept after assembly (default: ${MIN_LENGTH})
  --tmp-dir DIR Temp directory for Bakta (default: ${TMP_DIR})
  -g SIZE       Genome size for Flye (long_np/long_pb only, e.g. 5.8m).
                Only takes effect together with --asm-coverage; unset by
                default (Flye uses all reads for the initial assembly).
  --asm-coverage N
                Cap Flye's initial disjointig-assembly coverage to N x the
                longest reads (long_np/long_pb only; requires -g). Fixes
                Flye's "No disjointigs were assembled" failure on very
                high-coverage runs (documented at >800x, seen in practice
                around 600x+); 40-50 is Flye's own recommended range.
  --skip-checkm Skip the CheckM step
  --skip-bakta  Skip the Bakta annotation step
  --force       Re-run steps even if their output already exists
  -h            Show this help and exit

short:     Trimmomatic -> SPAdes -> reformat.sh -> QUAST -> BBMap
long_np:   Filtlong -> Flye --nano-* -> medaka -> reformat.sh -> QUAST -> minimap2/mosdepth
long_pb:   Filtlong -> Flye --pacbio-hifi -> reformat.sh -> QUAST -> minimap2/mosdepth
hybrid_np: Trimmomatic + Filtlong -> SPAdes --nanopore     -> reformat.sh -> QUAST -> BBMap + minimap2/mosdepth (separate)
hybrid_pb: Trimmomatic + Filtlong -> SPAdes -s (HiFi as extra library) -> reformat.sh -> QUAST -> BBMap + minimap2/mosdepth (separate)

Once every sample has an assembly, CheckM runs once over all of them, then
Bakta annotates each sample's final assembly (same for every mode).

Long/hybrid config (see config.sh): FILTLONG_MIN_LENGTH, FILTLONG_KEEP_PERCENT,
FLYE_READ_TYPE (long_np), FLYE_PACBIO_READ_TYPE (long_pb), MEDAKA_MODEL (long_np only).
EOF
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    -r) RAW_DIR="$2"; shift 2 ;;
    -l) LONG_DIR="$2"; shift 2 ;;
    -o) OUT_DIR="$2"; shift 2 ;;
    -s) SAMPLES="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -a) ADAPTER="$2"; shift 2 ;;
    -d) BAKTA_DB="$2"; shift 2 ;;
    -m) MIN_LENGTH="$2"; shift 2 ;;
    --tmp-dir) TMP_DIR="$2"; shift 2 ;;
    -g) GENOME_SIZE="$2"; shift 2 ;;
    --asm-coverage) ASM_COVERAGE="$2"; shift 2 ;;
    --skip-checkm) SKIP_CHECKM=1; shift ;;
    --skip-bakta) SKIP_BAKTA=1; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$MODE" ]]; then
  echo "ERROR: --mode is required (short, long_np, long_pb, hybrid_np, hybrid_pb)" >&2
  usage
  exit 1
fi
case "$MODE" in
  short|long_np|long_pb|hybrid_np|hybrid_pb) ;;
  *) echo "Invalid --mode '$MODE' (expected short, long_np, long_pb, hybrid_np, or hybrid_pb)" >&2; exit 1 ;;
esac
USE_SHORT=0; USE_LONG=0; PLATFORM=""
[[ "$MODE" == "short" || "$MODE" == hybrid_* ]] && USE_SHORT=1
[[ "$MODE" == long_* || "$MODE" == hybrid_* ]] && USE_LONG=1
case "$MODE" in
  long_np|hybrid_np) PLATFORM="nanopore"; MINIMAP_PRESET="map-ont" ;;
  long_pb|hybrid_pb) PLATFORM="pacbio"; MINIMAP_PRESET="map-hifi" ;;
esac

# ---------------------------------------------------------------------------
# Derived paths (numbered so `ls` sorts them in pipeline order)
# ---------------------------------------------------------------------------
LOG_DIR="$OUT_DIR/00_logs"
QC_DIR="$OUT_DIR/01_QC"
ASSEMBLY_DIR="$OUT_DIR/02_assembly"
QUAST_DIR="$OUT_DIR/03_quast"
CHECKM_BIN_DIR="$OUT_DIR/04_checkm/bins"
CHECKM_OUT_DIR="$OUT_DIR/04_checkm/output"
BAKTA_DIR="$OUT_DIR/05_bakta"

mkdir -p "$QC_DIR" "$ASSEMBLY_DIR" "$QUAST_DIR" "$CHECKM_BIN_DIR" "$BAKTA_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
if [[ -t 1 ]]; then
  C_RESET=$'\033[0m'; C_BOLD=$'\033[1m'; C_DIM=$'\033[2m'
  C_RED=$'\033[31m'; C_GREEN=$'\033[32m'; C_YELLOW=$'\033[33m'
  C_BLUE=$'\033[34m'; C_CYAN=$'\033[36m'
else
  C_RESET=""; C_BOLD=""; C_DIM=""
  C_RED=""; C_GREEN=""; C_YELLOW=""; C_BLUE=""; C_CYAN=""
fi

log() {
  echo "${C_DIM}[$(date '+%Y-%m-%d %H:%M:%S')]${C_RESET} $*"
}

banner() {
  local title=" $* "
  local line
  line=$(printf '%*s' "${#title}" '' | tr ' ' '=')
  echo "${C_BOLD}${C_BLUE}${line}${C_RESET}"
  echo "${C_BOLD}${C_BLUE}${title}${C_RESET}"
  echo "${C_BOLD}${C_BLUE}${line}${C_RESET}"
}

format_duration() {
  local s=$1
  printf '%02dh%02dm%02ds' $((s / 3600)) $((s % 3600 / 60)) $((s % 60))
}

check_deps() {
  local tools=(reformat.sh quast.py)
  [[ "$USE_SHORT" -eq 1 ]] && tools+=(trimmomatic spades.py bbmap.sh)
  [[ "$USE_LONG" -eq 1 ]] && tools+=(filtlong minimap2 samtools mosdepth)
  [[ "$MODE" == long_* ]] && tools+=(flye)
  [[ "$MODE" == "long_np" ]] && tools+=(medaka_consensus)
  [[ "$SKIP_CHECKM" -eq 1 ]] || tools+=(checkm)
  [[ "$SKIP_BAKTA" -eq 1 ]] || tools+=(bakta)
  local missing=()
  for t in "${tools[@]}"; do
    command -v "$t" >/dev/null 2>&1 || missing+=("$t")
  done
  if [[ ${#missing[@]} -gt 0 ]]; then
    echo "Missing required tool(s) on PATH: ${missing[*]}" >&2
    exit 1
  fi
}

# Catch config path typos (e.g. a missing leading '/') before any sample
# runs, instead of failing on the last pipeline step after hours of work.
check_config() {
  local errors=()
  if [[ "$USE_SHORT" -eq 1 ]]; then
    [[ "$ADAPTER" == /* ]] || errors+=("ADAPTER is not an absolute path: '$ADAPTER'")
    [[ -f "$ADAPTER" ]] || errors+=("ADAPTER file not found: '$ADAPTER'")
  fi
  if [[ "$MODE" == "long_np" ]]; then
    [[ -n "$MEDAKA_MODEL" ]] || errors+=("MEDAKA_MODEL is not set (see config.sh)")
  fi
  if [[ -n "$ASM_COVERAGE" && -z "$GENOME_SIZE" ]]; then
    errors+=("--asm-coverage requires -g (genome size) to also be set")
  fi
  if [[ -n "$GENOME_SIZE" && -z "$ASM_COVERAGE" ]]; then
    errors+=("-g (genome size) has no effect without --asm-coverage")
  fi
  if [[ "$SKIP_BAKTA" -eq 0 ]]; then
    [[ "$BAKTA_DB" == /* ]] || errors+=("BAKTA_DB is not an absolute path: '$BAKTA_DB'")
    [[ -d "$BAKTA_DB" ]] || errors+=("BAKTA_DB directory not found: '$BAKTA_DB'")
    [[ "$TMP_DIR" == /* ]] || errors+=("TMP_DIR is not an absolute path: '$TMP_DIR'")
  fi
  if [[ ${#errors[@]} -gt 0 ]]; then
    printf 'Config error: %s\n' "${errors[@]}" >&2
    echo "Fix these in config.sh (or pass the equivalent flags) before running." >&2
    exit 1
  fi
}

# Run a step's command, logging to a per-sample/per-step file.
# Skips the command if $3 (marker/output file) already exists, unless --force.
STEP_RESULTS=()

run_step() {
  local name="$1" logfile="$2" marker="$3"
  shift 3
  if [[ "$FORCE" -eq 0 && -n "$marker" && -e "$marker" ]]; then
    log "${C_YELLOW}${C_BOLD}SKIP ${C_RESET} ${name} ${C_DIM}(already exists: $marker)${C_RESET}"
    STEP_RESULTS+=("SKIP  $name")
    return 0
  fi
  log "${C_CYAN}${C_BOLD}START${C_RESET} ${name}"
  local start_ts elapsed
  start_ts=$(date +%s)
  if ! "$@" >"$logfile" 2>&1; then
    log "${C_RED}${C_BOLD}FAIL ${C_RESET} ${name} ${C_DIM}(see $logfile)${C_RESET}"
    STEP_RESULTS+=("FAIL  $name")
    return 1
  fi
  elapsed=$(( $(date +%s) - start_ts ))
  log "${C_GREEN}${C_BOLD}DONE ${C_RESET} ${name} ${C_DIM}($(format_duration "$elapsed"))${C_RESET}"
  STEP_RESULTS+=("DONE  $name")
}

# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------
# Short-read pairs are recognized as <sample>_1/_2.fastq.gz,
# <sample>_R1/_R2.fastq.gz, or Illumina's own bcl2fastq/bcl-convert default
# <sample>_S#_L#_R1_001.fastq.gz/_R2_001.fastq.gz (all common). Sets
# RAW_R1/RAW_R2, or leaves them empty if no convention matches.
resolve_short_reads() {
  local sample="$1"
  if [[ -f "$RAW_DIR/${sample}_1.fastq.gz" ]]; then
    RAW_R1="$RAW_DIR/${sample}_1.fastq.gz"
    RAW_R2="$RAW_DIR/${sample}_2.fastq.gz"
  elif [[ -f "$RAW_DIR/${sample}_R1.fastq.gz" ]]; then
    RAW_R1="$RAW_DIR/${sample}_R1.fastq.gz"
    RAW_R2="$RAW_DIR/${sample}_R2.fastq.gz"
  else
    RAW_R1=""; RAW_R2=""
    local r1_matches=("$RAW_DIR/${sample}"_S*_L*_R1_001.fastq.gz)
    if [[ -e "${r1_matches[0]}" ]]; then
      if [[ ${#r1_matches[@]} -gt 1 ]]; then
        echo "ERROR: sample '$sample' has multiple lane files matching" \
             "${sample}_S*_L*_R1_001.fastq.gz: ${r1_matches[*]}." \
             "This script doesn't merge multi-lane runs automatically --" \
             "concatenate each read direction's lanes into one file first, e.g.:" \
             "  cat ${sample}_S*_L*_R1_001.fastq.gz > ${sample}_1.fastq.gz" \
             "  cat ${sample}_S*_L*_R2_001.fastq.gz > ${sample}_2.fastq.gz" \
             "(then re-run; the _1/_2 pair above takes priority over this pattern)." >&2
        exit 1
      fi
      RAW_R1="${r1_matches[0]}"
      RAW_R2="${RAW_R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    fi
  fi
}

# Explicit =() (not just `declare -a`): an array that's declared but never
# assigned an element still trips `set -u` on `${arr[@]}`/`${#arr[@]}` in
# bash, even in recent versions.
SAMPLE_LIST=()
if [[ -n "$SAMPLES" ]]; then
  IFS=',' read -r -a SAMPLE_LIST <<< "$SAMPLES"
elif [[ "$USE_SHORT" -eq 1 ]]; then
  for f in "$RAW_DIR"/*_1.fastq.gz "$RAW_DIR"/*_R1.fastq.gz "$RAW_DIR"/*_S*_L*_R1_001.fastq.gz; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    base="${base%_1.fastq.gz}"
    base="${base%_R1.fastq.gz}"
    base="${base%_S*_L*_R1_001.fastq.gz}"
    SAMPLE_LIST+=("$base")
  done
  # Dedupe in case a sample somehow matched more than one pattern, or has
  # reads split across multiple lanes (resolve_short_reads() errors out on
  # that case below rather than silently dropping lanes).
  if [[ ${#SAMPLE_LIST[@]} -gt 0 ]]; then
    mapfile -t SAMPLE_LIST < <(printf '%s\n' "${SAMPLE_LIST[@]}" | sort -u)
  fi
else
  for f in "$LONG_DIR"/*.fastq.gz; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    SAMPLE_LIST+=("${base%.fastq.gz}")
  done
fi

if [[ ${#SAMPLE_LIST[@]} -eq 0 ]]; then
  echo "No samples found (mode=$MODE; checked '$RAW_DIR' / '$LONG_DIR')" >&2
  exit 1
fi

for sample in "${SAMPLE_LIST[@]}"; do
  if [[ "$USE_SHORT" -eq 1 ]]; then
    resolve_short_reads "$sample"
    if [[ -z "$RAW_R1" || ! -f "$RAW_R1" || ! -f "$RAW_R2" ]]; then
      echo "Missing short reads for sample '$sample' (looked for" \
           "${sample}_1.fastq.gz/_2.fastq.gz and ${sample}_R1.fastq.gz/_R2.fastq.gz" \
           "in '$RAW_DIR')" >&2
      exit 1
    fi
  fi
  if [[ "$USE_LONG" -eq 1 ]]; then
    lr="$LONG_DIR/${sample}.fastq.gz"
    if [[ ! -f "$lr" ]]; then
      echo "Missing long reads for sample '$sample' ($lr)" >&2
      exit 1
    fi
  fi
done

check_deps
check_config

# Only takes effect for long_np/long_pb; empty unless both -g and
# --asm-coverage were passed (validated together in check_config).
FLYE_EXTRA_ARGS=()
if [[ -n "$GENOME_SIZE" && -n "$ASM_COVERAGE" ]]; then
  FLYE_EXTRA_ARGS=(--genome-size "$GENOME_SIZE" --asm-coverage "$ASM_COVERAGE")
fi

banner "Bacteria WGS Pipeline"
echo "  ${C_BOLD}Mode${C_RESET}      : $MODE${PLATFORM:+ ($PLATFORM)}"
echo "  ${C_BOLD}Samples${C_RESET}   : ${SAMPLE_LIST[*]}"
[[ "$USE_SHORT" -eq 1 ]] && echo "  ${C_BOLD}Short reads${C_RESET}: $RAW_DIR"
[[ "$USE_LONG" -eq 1 ]] && echo "  ${C_BOLD}Long reads${C_RESET} : $LONG_DIR"
echo "  ${C_BOLD}Output${C_RESET}    : $OUT_DIR"
echo "  ${C_BOLD}Threads${C_RESET}   : $THREADS"
echo

# ---------------------------------------------------------------------------
# Per-sample steps: QC/filtering -> assembly -> assembly QC
# ---------------------------------------------------------------------------
for sample in "${SAMPLE_LIST[@]}"; do
  echo
  echo "${C_BOLD}${C_BLUE}--- Sample: ${sample} (mode: $MODE) ---${C_RESET}"
  sample_qc_dir="$QC_DIR/$sample"
  sample_asm_dir="$ASSEMBLY_DIR/$sample"
  sample_quast_dir="$QUAST_DIR/$sample"
  mkdir -p "$sample_qc_dir" "$sample_quast_dir"

  # --- QC / filtering ---
  if [[ "$USE_SHORT" -eq 1 ]]; then
    resolve_short_reads "$sample"
    filt_r1="$sample_qc_dir/${sample}_filtered_1.fastq.gz"
    filt_r2="$sample_qc_dir/${sample}_filtered_2.fastq.gz"
    run_step "trimmomatic:$sample" "$LOG_DIR/${sample}_trimmomatic.log" "$filt_r1" \
      trimmomatic PE -threads "$THREADS" -phred33 \
        -summary "$sample_qc_dir/summary.txt" -compressLevel 5 \
        "$RAW_R1" "$RAW_R2" \
        "$filt_r1" "$sample_qc_dir/${sample}_unpaired_1.fastq.gz" \
        "$filt_r2" "$sample_qc_dir/${sample}_unpaired_2.fastq.gz" \
        ILLUMINACLIP:"${ADAPTER}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  fi

  if [[ "$USE_LONG" -eq 1 ]]; then
    filt_long="$sample_qc_dir/${sample}_long_filtered.fastq.gz"
    run_step "filtlong:$sample" "$LOG_DIR/${sample}_filtlong.log" "$filt_long" \
      bash -c 'set -o pipefail; filtlong --min_length "$1" --keep_percent "$2" "$3" | gzip > "$4"' _ \
        "$FILTLONG_MIN_LENGTH" "$FILTLONG_KEEP_PERCENT" "$LONG_DIR/${sample}.fastq.gz" "$filt_long"
  fi

  # --- Assembly ---
  case "$MODE" in
    short)
      run_step "spades:$sample" "$LOG_DIR/${sample}_spades.log" "$sample_asm_dir/scaffolds.fasta" \
        spades.py -1 "$filt_r1" -2 "$filt_r2" -o "$sample_asm_dir" \
          --isolate --cov-cutoff auto -k 21,33,55,77,99,127 -t "$THREADS"
      raw_assembly="$sample_asm_dir/scaffolds.fasta"
      ;;
    hybrid_np)
      run_step "spades_hybrid:$sample" "$LOG_DIR/${sample}_spades.log" "$sample_asm_dir/scaffolds.fasta" \
        spades.py -1 "$filt_r1" -2 "$filt_r2" --nanopore "$filt_long" -o "$sample_asm_dir" \
          --isolate --cov-cutoff auto -k 21,33,55,77,99,127 -t "$THREADS"
      raw_assembly="$sample_asm_dir/scaffolds.fasta"
      ;;
    hybrid_pb)
      # SPAdes has no --pacbio-hifi flag: --pacbio is for noisy PacBio CLR
      # reads. HiFi/CCS reads are accurate enough to feed in as an extra
      # single-read library via -s instead.
      run_step "spades_hybrid:$sample" "$LOG_DIR/${sample}_spades.log" "$sample_asm_dir/scaffolds.fasta" \
        spades.py -1 "$filt_r1" -2 "$filt_r2" -s "$filt_long" -o "$sample_asm_dir" \
          --isolate --cov-cutoff auto -k 21,33,55,77,99,127 -t "$THREADS"
      raw_assembly="$sample_asm_dir/scaffolds.fasta"
      ;;
    long_np)
      run_step "flye:$sample" "$LOG_DIR/${sample}_flye.log" "$sample_asm_dir/assembly.fasta" \
        flye "$FLYE_READ_TYPE" "$filt_long" --out-dir "$sample_asm_dir" --threads "$THREADS" \
          "${FLYE_EXTRA_ARGS[@]}"

      medaka_dir="$sample_asm_dir/medaka"
      run_step "medaka:$sample" "$LOG_DIR/${sample}_medaka.log" "$medaka_dir/consensus.fasta" \
        medaka_consensus -i "$filt_long" -d "$sample_asm_dir/assembly.fasta" \
          -o "$medaka_dir" -t "$THREADS" -m "$MEDAKA_MODEL"
      raw_assembly="$medaka_dir/consensus.fasta"
      ;;
    long_pb)
      run_step "flye:$sample" "$LOG_DIR/${sample}_flye.log" "$sample_asm_dir/assembly.fasta" \
        flye "$FLYE_PACBIO_READ_TYPE" "$filt_long" --out-dir "$sample_asm_dir" --threads "$THREADS" \
          "${FLYE_EXTRA_ARGS[@]}"
      # No medaka: it's an ONT-trained polisher, not applicable to PacBio HiFi.
      raw_assembly="$sample_asm_dir/assembly.fasta"
      ;;
  esac

  final_fasta="$sample_asm_dir/scaffolds_final.fasta"
  run_step "reformat:$sample" "$LOG_DIR/${sample}_reformat.log" "$final_fasta" \
    reformat.sh in="$raw_assembly" out="$final_fasta" minlength="$MIN_LENGTH"

  # QUAST gets the relevant reads per mode so its report includes
  # mapping-based metrics (mapped %, coverage, mismatches/indels), not
  # just bare assembly stats.
  quast_read_args=()
  [[ "$USE_SHORT" -eq 1 ]] && quast_read_args+=(--pe1 "$filt_r1" --pe2 "$filt_r2")
  if [[ "$USE_LONG" -eq 1 ]]; then
    if [[ "$PLATFORM" == "pacbio" ]]; then
      quast_read_args+=(--pacbio "$filt_long")
    else
      quast_read_args+=(--nanopore "$filt_long")
    fi
  fi
  run_step "quast:$sample" "$LOG_DIR/${sample}_quast.log" "$sample_quast_dir/report.txt" \
    quast.py "$final_fasta" -o "$sample_quast_dir" -t "$THREADS" "${quast_read_args[@]}"

  # --- Coverage ---
  case "$MODE" in
    short)
      run_step "bbmap_coverage:$sample" "$LOG_DIR/${sample}_coverage.log" "$sample_quast_dir/coverage_stats.txt" \
        bbmap.sh ref="$final_fasta" in1="$filt_r1" in2="$filt_r2" \
          covstats="$sample_quast_dir/coverage_stats.txt" threads="$THREADS"
      ;;
    long_np|long_pb)
      run_step "coverage:$sample" "$LOG_DIR/${sample}_coverage.log" \
          "$sample_quast_dir/${sample}_coverage.mosdepth.summary.txt" \
        env THREADS="$THREADS" FASTA="$final_fasta" LONG="$filt_long" PRESET="$MINIMAP_PRESET" \
            BAM="$sample_quast_dir/${sample}_long.sorted.bam" \
            PREFIX="$sample_quast_dir/${sample}_coverage" \
        bash -c '
          set -euo pipefail
          minimap2 -ax "$PRESET" -t "$THREADS" "$FASTA" "$LONG" \
            | samtools sort -@ "$THREADS" -o "$BAM" -
          samtools index "$BAM"
          mosdepth -n -t "$THREADS" --fast-mode "$PREFIX" "$BAM"
        '
      ;;
    hybrid_np|hybrid_pb)
      # Per-platform numbers (always useful on their own), plus one true
      # combined figure via samtools merge. The earlier merge failure was
      # most likely because the short-read bam was never coordinate-sorted
      # while the long-read bam was -- samtools merge requires matching
      # sort order across inputs. Sorting both the same way fixes that.
      # Marker is the sorted+indexed bam (what total_coverage below needs),
      # not the covstats file -- an older pipeline version only produced
      # the covstats file here, so using it as the marker would let a
      # resumed run skip this step while still missing the bam.
      run_step "bbmap_coverage:$sample" "$LOG_DIR/${sample}_coverage.log" \
          "$sample_quast_dir/${sample}_illumina.sorted.bam.bai" \
        env THREADS="$THREADS" FASTA="$final_fasta" R1="$filt_r1" R2="$filt_r2" \
            RAW_BAM="$sample_quast_dir/${sample}_illumina.unsorted.bam" \
            BAM="$sample_quast_dir/${sample}_illumina.sorted.bam" \
            COV="$sample_quast_dir/${sample}_illumina_cov.txt" \
        bash -c '
          set -euo pipefail
          bbmap.sh ref="$FASTA" in1="$R1" in2="$R2" out="$RAW_BAM" covstats="$COV" threads="$THREADS"
          samtools sort -@ "$THREADS" -o "$BAM" "$RAW_BAM"
          samtools index "$BAM"
          rm -f "$RAW_BAM"
        '

      run_step "long_coverage:$sample" "$LOG_DIR/${sample}_long_coverage.log" \
          "$sample_quast_dir/${sample}_long_coverage.mosdepth.summary.txt" \
        env THREADS="$THREADS" FASTA="$final_fasta" LONG="$filt_long" PRESET="$MINIMAP_PRESET" \
            BAM="$sample_quast_dir/${sample}_long.sorted.bam" \
            PREFIX="$sample_quast_dir/${sample}_long_coverage" \
        bash -c '
          set -euo pipefail
          minimap2 -ax "$PRESET" -t "$THREADS" "$FASTA" "$LONG" \
            | samtools sort -@ "$THREADS" -o "$BAM" -
          samtools index "$BAM"
          mosdepth -n -t "$THREADS" --fast-mode "$PREFIX" "$BAM"
        '

      run_step "total_coverage:$sample" "$LOG_DIR/${sample}_total_coverage.log" \
          "$sample_quast_dir/${sample}_total_coverage.mosdepth.summary.txt" \
        env THREADS="$THREADS" \
            ILL_BAM="$sample_quast_dir/${sample}_illumina.sorted.bam" \
            LONG_BAM="$sample_quast_dir/${sample}_long.sorted.bam" \
            MERGED_BAM="$sample_quast_dir/${sample}_merged.sorted.bam" \
            PREFIX="$sample_quast_dir/${sample}_total_coverage" \
        bash -c '
          set -euo pipefail
          samtools merge -f -@ "$THREADS" "$MERGED_BAM" "$ILL_BAM" "$LONG_BAM"
          samtools index "$MERGED_BAM"
          mosdepth -n -t "$THREADS" --fast-mode "$PREFIX" "$MERGED_BAM"
        '
      ;;
  esac

  cp -f "$final_fasta" "$CHECKM_BIN_DIR/${sample}_final.fasta"
done

# ---------------------------------------------------------------------------
# CheckM (batch, over every sample's final assembly at once)
# ---------------------------------------------------------------------------
if [[ "$SKIP_CHECKM" -eq 0 ]]; then
  run_step "checkm_lineage_wf" "$LOG_DIR/checkm_lineage_wf.log" "$CHECKM_OUT_DIR/lineage.ms" \
    checkm lineage_wf "$CHECKM_BIN_DIR" "$CHECKM_OUT_DIR" -x fasta --tab_table -t "$THREADS"

  run_step "checkm_qa" "$LOG_DIR/checkm_qa.log" "$CHECKM_OUT_DIR/checkm_results.txt" \
    checkm qa -o 2 "$CHECKM_OUT_DIR/lineage.ms" "$CHECKM_OUT_DIR" \
      --tab_table -f "$CHECKM_OUT_DIR/checkm_results.txt"
else
  log "${C_YELLOW}${C_BOLD}SKIP ${C_RESET} CheckM ${C_DIM}(--skip-checkm)${C_RESET}"
  STEP_RESULTS+=("SKIP  CheckM")
fi

# ---------------------------------------------------------------------------
# Bakta annotation (per sample, from the CheckM bin copy of the assembly)
# ---------------------------------------------------------------------------
if [[ "$SKIP_BAKTA" -eq 0 ]]; then
  mkdir -p "$TMP_DIR"
  for sample in "${SAMPLE_LIST[@]}"; do
    run_step "bakta:$sample" "$LOG_DIR/${sample}_bakta.log" "$BAKTA_DIR/$sample/${sample}_final.gff3" \
      bakta --db "$BAKTA_DB" -t "$THREADS" --force \
        --output "$BAKTA_DIR/$sample" "$CHECKM_BIN_DIR/${sample}_final.fasta" \
        --tmp-dir "$TMP_DIR"
  done
else
  log "${C_YELLOW}${C_BOLD}SKIP ${C_RESET} Bakta ${C_DIM}(--skip-bakta)${C_RESET}"
  STEP_RESULTS+=("SKIP  Bakta")
fi

echo
banner "Summary"
for r in "${STEP_RESULTS[@]}"; do
  status="${r%% *}"
  case "$status" in
    DONE) color="$C_GREEN" ;;
    SKIP) color="$C_YELLOW" ;;
    FAIL) color="$C_RED" ;;
    *) color="$C_RESET" ;;
  esac
  echo "  ${color}${C_BOLD}${r}${C_RESET}"
done
echo
log "${C_GREEN}${C_BOLD}Pipeline finished.${C_RESET} Results under: $OUT_DIR"
