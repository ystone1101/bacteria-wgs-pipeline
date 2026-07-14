// short mode: Illumina reads only.
process QUAST_SHORT {
    tag "$sample"
    publishDir "${params.outdir}/03_quast/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(r1), path(r2)

    output:
    tuple val(sample), path("report.txt"), emit: report
    path "*"

    script:
    """
    quast.py ${assembly} -o . -t ${task.cpus} --pe1 ${r1} --pe2 ${r2}
    """
}
