#!/usr/bin/env nextflow
//
// Personal Nextflow DSL2 port of bacteria_wgs_pipeline.sh
// Modes (params.mode, required, no default): short, long_np, long_pb,
// hybrid_np, hybrid_pb -- same tool chains as the bash script, see
// nextflow.config / README for details.
//
// Usage:
//   nextflow run main.nf -profile conda --mode short     --input samplesheet.csv --adapter /path/adapters.fa --bakta_db /path/db
//   nextflow run main.nf -profile conda --mode long_np    --input samplesheet.csv --medaka_model r941_prom_hac_g507 --bakta_db /path/db
//   nextflow run main.nf -profile conda --mode hybrid_pb  --input samplesheet.csv --adapter /path/adapters.fa --bakta_db /path/db
//
// samplesheet.csv columns: sample,fastq_1,fastq_2,long_reads
// (leave fastq_1/fastq_2 blank for long_* modes, long_reads blank for short)

nextflow.enable.dsl = 2

include { TRIMMOMATIC }           from './modules/local/trimmomatic'
include { FILTLONG }              from './modules/local/filtlong'
include { SPADES_SHORT }          from './modules/local/spades_short'
include { SPADES_HYBRID }         from './modules/local/spades_hybrid'
include { FLYE }                  from './modules/local/flye'
include { MEDAKA }                from './modules/local/medaka'
include { REFORMAT }              from './modules/local/reformat'
include { QUAST_SHORT }           from './modules/local/quast_short'
include { QUAST_LONG }            from './modules/local/quast_long'
include { QUAST_HYBRID }          from './modules/local/quast_hybrid'
include { COVERAGE_SHORT }        from './modules/local/coverage_short'
include { COVERAGE_SHORT_HYBRID } from './modules/local/coverage_short_hybrid'
include { COVERAGE_LONG }         from './modules/local/coverage_long'
include { COVERAGE_TOTAL }        from './modules/local/coverage_total'
include { CHECKM }                from './modules/local/checkm'
include { BAKTA }                 from './modules/local/bakta'

def VALID_MODES = ['short', 'long_np', 'long_pb', 'hybrid_np', 'hybrid_pb']

workflow {

  // --- Preflight validation (mirrors check_config() in the bash script) ---
  if (!params.mode) {
    error "ERROR: --mode is required (${VALID_MODES.join(', ')})"
  }
  if (!(params.mode in VALID_MODES)) {
    error "ERROR: invalid --mode '${params.mode}' (expected one of: ${VALID_MODES.join(', ')})"
  }
  if (!params.input) {
    error "ERROR: --input <samplesheet.csv> is required"
  }

  def use_short = (params.mode == 'short' || params.mode.startsWith('hybrid'))
  def use_long  = (params.mode.startsWith('long') || params.mode.startsWith('hybrid'))
  def platform  = params.mode.endsWith('_np') ? 'nanopore' : (params.mode.endsWith('_pb') ? 'pacbio' : null)
  def minimap_preset     = (platform == 'pacbio') ? 'map-hifi' : 'map-ont'
  def quast_platform_flag = (platform == 'pacbio') ? '--pacbio' : '--nanopore'

  if (use_short && !params.adapter) {
    error "ERROR: --adapter is required for mode ${params.mode}"
  }
  if (params.mode == 'long_np' && !params.medaka_model) {
    error "ERROR: --medaka_model is required for mode long_np (see README for picking one)"
  }
  if (!params.skip_bakta && !params.bakta_db) {
    error "ERROR: --bakta_db is required unless --skip_bakta is set"
  }

  // --- Samplesheet ---
  ch_rows = Channel.fromPath(params.input, checkIfExists: true)
    .splitCsv(header: true)

  if (use_short) {
    ch_short_raw = ch_rows.map { row ->
      if (!row.fastq_1 || !row.fastq_2) {
        error "ERROR: sample '${row.sample}' is missing fastq_1/fastq_2 (required for mode ${params.mode})"
      }
      tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
    }
  }
  if (use_long) {
    ch_long_raw = ch_rows.map { row ->
      if (!row.long_reads) {
        error "ERROR: sample '${row.sample}' is missing long_reads (required for mode ${params.mode})"
      }
      tuple(row.sample, file(row.long_reads, checkIfExists: true))
    }
  }

  // --- QC / filtering ---
  ch_short_filtered = use_short ? TRIMMOMATIC(ch_short_raw, file(params.adapter, checkIfExists: true)).reads : Channel.empty()
  ch_long_filtered  = use_long ? FILTLONG(ch_long_raw).reads : Channel.empty()

  // --- Assembly (per mode) ---
  switch (params.mode) {
    case 'short':
      SPADES_SHORT(ch_short_filtered)
      ch_raw_assembly = SPADES_SHORT.out.assembly
      break
    case 'hybrid_np':
      SPADES_HYBRID(ch_short_filtered.join(ch_long_filtered), '--nanopore')
      ch_raw_assembly = SPADES_HYBRID.out.assembly
      break
    case 'hybrid_pb':
      // SPAdes has no --pacbio-hifi flag; HiFi reads go in as an extra
      // single-read library via -s (see bash script's comment).
      SPADES_HYBRID(ch_short_filtered.join(ch_long_filtered), '-s')
      ch_raw_assembly = SPADES_HYBRID.out.assembly
      break
    case 'long_np':
      FLYE(ch_long_filtered, params.flye_read_type)
      MEDAKA(FLYE.out.assembly.join(ch_long_filtered))
      ch_raw_assembly = MEDAKA.out.assembly
      break
    case 'long_pb':
      // No medaka: it's an ONT-trained polisher, not applicable to HiFi.
      FLYE(ch_long_filtered, params.flye_pacbio_read_type)
      ch_raw_assembly = FLYE.out.assembly
      break
  }

  REFORMAT(ch_raw_assembly)
  ch_final_assembly = REFORMAT.out.assembly

  // --- QUAST ---
  if (params.mode == 'short') {
    QUAST_SHORT(ch_final_assembly.join(ch_short_filtered))
  } else if (use_short && use_long) {
    QUAST_HYBRID(ch_final_assembly.join(ch_short_filtered).join(ch_long_filtered), quast_platform_flag)
  } else {
    QUAST_LONG(ch_final_assembly.join(ch_long_filtered), quast_platform_flag)
  }

  // --- Coverage ---
  if (params.mode == 'short') {
    COVERAGE_SHORT(ch_final_assembly.join(ch_short_filtered))
  } else if (use_short && use_long) {
    COVERAGE_SHORT_HYBRID(ch_final_assembly.join(ch_short_filtered))
    COVERAGE_LONG(ch_final_assembly.join(ch_long_filtered), minimap_preset, 'long_coverage')
    COVERAGE_TOTAL(COVERAGE_SHORT_HYBRID.out.bam.join(COVERAGE_LONG.out.bam))
  } else {
    COVERAGE_LONG(ch_final_assembly.join(ch_long_filtered), minimap_preset, 'coverage')
  }

  // --- CheckM (batch, over every sample's final assembly at once) ---
  if (!params.skip_checkm) {
    CHECKM(ch_final_assembly.map { it[1] }.collect())
  }

  // --- Bakta annotation (per sample) ---
  if (!params.skip_bakta) {
    BAKTA(ch_final_assembly)
  }
}
