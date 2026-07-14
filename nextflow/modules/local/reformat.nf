// Output is named "${sample}_final.fasta" (not a generic name) so that
// CHECKM can .collect() every sample's assembly into one bin directory
// without filename collisions.
process REFORMAT {
    tag "$sample"
    publishDir "${params.outdir}/02_assembly/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(raw_assembly)

    output:
    tuple val(sample), path("${sample}_final.fasta"), emit: assembly

    script:
    """
    reformat.sh in=${raw_assembly} out=${sample}_final.fasta minlength=${params.min_length}
    """
}
