// Batch step: runs once over every sample's final assembly at once, same
// as the bash pipeline. Assumes `checkm data setRoot <dir>` has already
// been run once for this user (see README) -- that's a one-time global
// config, not something this process can set up itself.
process CHECKM {
    publishDir "${params.outdir}/04_checkm/output", mode: 'copy'

    input:
    path bins

    output:
    path "checkm_results.txt", emit: results
    path "*"

    script:
    """
    mkdir -p bins
    mv ${bins} bins/
    checkm lineage_wf bins . -x fasta --tab_table -t ${task.cpus}
    checkm qa -o 2 lineage.ms . --tab_table -f checkm_results.txt
    """
}
