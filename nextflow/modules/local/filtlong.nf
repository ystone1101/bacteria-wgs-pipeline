process FILTLONG {
    tag "$sample"
    publishDir "${params.outdir}/01_QC/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(long_reads)

    output:
    tuple val(sample), path("${sample}_long_filtered.fastq.gz"), emit: reads

    script:
    """
    filtlong --min_length ${params.filtlong_min_length} \\
      --keep_percent ${params.filtlong_keep_percent} \\
      ${long_reads} | gzip > ${sample}_long_filtered.fastq.gz
    """
}
