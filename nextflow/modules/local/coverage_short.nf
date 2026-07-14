// short mode only: BBMap covstats against the final assembly.
process COVERAGE_SHORT {
    tag "$sample"
    publishDir "${params.outdir}/03_quast/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(r1), path(r2)

    output:
    tuple val(sample), path("coverage_stats.txt"), emit: stats

    script:
    """
    bbmap.sh ref=${assembly} in1=${r1} in2=${r2} covstats=coverage_stats.txt threads=${task.cpus}
    """
}
