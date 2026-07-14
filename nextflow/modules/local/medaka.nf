// long_np only. medaka is ONT-trained -- not run for long_pb (HiFi).
process MEDAKA {
    tag "$sample"
    publishDir "${params.outdir}/02_assembly/${sample}/medaka", mode: 'copy'

    input:
    tuple val(sample), path(draft), path(long_reads)

    output:
    tuple val(sample), path("consensus.fasta"), emit: assembly

    script:
    """
    medaka_consensus -i ${long_reads} -d ${draft} -o . -t ${task.cpus} -m ${params.medaka_model}
    """
}
