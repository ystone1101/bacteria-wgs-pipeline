// hybrid_np / hybrid_pb: both short and long reads. platform_flag is
// '--nanopore' or '--pacbio'.
process QUAST_HYBRID {
    tag "$sample"
    publishDir "${params.outdir}/03_quast/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(r1), path(r2), path(long_reads)
    val platform_flag

    output:
    tuple val(sample), path("report.txt"), emit: report
    path "*"

    script:
    """
    quast.py ${assembly} -o . -t ${task.cpus} --pe1 ${r1} --pe2 ${r2} ${platform_flag} ${long_reads}
    """
}
