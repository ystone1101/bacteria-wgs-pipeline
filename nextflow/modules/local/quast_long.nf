// long_np / long_pb: long reads only. platform_flag is '--nanopore' or '--pacbio'.
process QUAST_LONG {
    tag "$sample"
    publishDir "${params.outdir}/03_quast/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(long_reads)
    val platform_flag

    output:
    tuple val(sample), path("report.txt"), emit: report
    path "*"

    script:
    """
    quast.py ${assembly} -o . -t ${task.cpus} ${platform_flag} ${long_reads}
    """
}
