// Shared by hybrid_np ('--nanopore') and hybrid_pb ('-s', SPAdes has no
// --pacbio-hifi flag -- HiFi reads go in as an extra single-read library).
process SPADES_HYBRID {
    tag "$sample"
    publishDir "${params.outdir}/02_assembly/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2), path(long_reads)
    val long_flag

    output:
    tuple val(sample), path("scaffolds.fasta"), emit: assembly
    path "spades.log", optional: true

    script:
    """
    spades.py -1 ${r1} -2 ${r2} ${long_flag} ${long_reads} -o . \\
      --isolate --cov-cutoff auto -k 21,33,55,77,99,127 -t ${task.cpus}
    """
}
