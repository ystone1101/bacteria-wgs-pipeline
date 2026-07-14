// hybrid_np / hybrid_pb only: merge the sorted Illumina and long-read bams
// into one combined-coverage figure.
process COVERAGE_TOTAL {
    tag "$sample"
    publishDir "${params.outdir}/03_quast/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(illumina_bam), path(illumina_bai), path(long_bam), path(long_bai)

    output:
    tuple val(sample), path("${sample}_total_coverage.mosdepth.summary.txt"), emit: summary
    path "${sample}_total_coverage.mosdepth*"
    path "${sample}_merged.sorted.bam*"

    script:
    """
    samtools merge -f -@ ${task.cpus} ${sample}_merged.sorted.bam ${illumina_bam} ${long_bam}
    samtools index ${sample}_merged.sorted.bam
    mosdepth -n -t ${task.cpus} --fast-mode ${sample}_total_coverage ${sample}_merged.sorted.bam
    """
}
