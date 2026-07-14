// long_np / long_pb (standalone) and hybrid_np / hybrid_pb (long side):
// minimap2 -> sorted/indexed bam -> mosdepth. preset is 'map-ont' or
// 'map-hifi'. prefix is 'coverage' for long_* modes, 'long_coverage' for
// hybrid_* modes (matches the bash script's two call sites); the bam
// filename itself is always "${sample}_long.sorted.bam" so hybrid's
// COVERAGE_TOTAL step can find it.
process COVERAGE_LONG {
    tag "$sample"
    publishDir "${params.outdir}/03_quast/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(long_reads)
    val preset
    val prefix

    output:
    tuple val(sample), path("${sample}_long.sorted.bam"), path("${sample}_long.sorted.bam.bai"), emit: bam
    tuple val(sample), path("${sample}_${prefix}.mosdepth.summary.txt"), emit: summary
    path "${sample}_${prefix}.mosdepth*"

    script:
    """
    minimap2 -ax ${preset} -t ${task.cpus} ${assembly} ${long_reads} \\
      | samtools sort -@ ${task.cpus} -o ${sample}_long.sorted.bam -
    samtools index ${sample}_long.sorted.bam
    mosdepth -n -t ${task.cpus} --fast-mode ${sample}_${prefix} ${sample}_long.sorted.bam
    """
}
