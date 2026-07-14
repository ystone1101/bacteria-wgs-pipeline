process FLYE {
    tag "$sample"
    publishDir "${params.outdir}/02_assembly/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(long_reads)
    val read_type_flag

    output:
    tuple val(sample), path("assembly.fasta"), emit: assembly
    path "flye.log", optional: true

    script:
    """
    flye ${read_type_flag} ${long_reads} --out-dir . --threads ${task.cpus}
    """
}
