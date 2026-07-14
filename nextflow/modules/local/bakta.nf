process BAKTA {
    tag "$sample"
    publishDir "${params.outdir}/05_bakta/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}_final.gff3"), emit: gff
    path "*"

    script:
    """
    mkdir -p tmp
    bakta --db ${params.bakta_db} -t ${task.cpus} --force \\
      --output . ${assembly} --tmp-dir tmp
    """
}
