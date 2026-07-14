// hybrid_np / hybrid_pb only: BBMap against the final assembly, then sort
// + index the resulting bam so it has the same sort order as the long-read
// bam from COVERAGE_LONG (samtools merge requires matching sort order).
process COVERAGE_SHORT_HYBRID {
    tag "$sample"
    publishDir "${params.outdir}/03_quast/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_illumina.sorted.bam"), path("${sample}_illumina.sorted.bam.bai"), emit: bam
    path "${sample}_illumina_cov.txt"

    script:
    """
    bbmap.sh ref=${assembly} in1=${r1} in2=${r2} out=${sample}_illumina.unsorted.bam \\
      covstats=${sample}_illumina_cov.txt threads=${task.cpus}
    samtools sort -@ ${task.cpus} -o ${sample}_illumina.sorted.bam ${sample}_illumina.unsorted.bam
    samtools index ${sample}_illumina.sorted.bam
    rm -f ${sample}_illumina.unsorted.bam
    """
}
