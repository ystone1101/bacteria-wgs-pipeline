process TRIMMOMATIC {
    tag "$sample"
    publishDir "${params.outdir}/01_QC/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)
    path adapter

    output:
    tuple val(sample), path("${sample}_filtered_1.fastq.gz"), path("${sample}_filtered_2.fastq.gz"), emit: reads
    path "${sample}_summary.txt"

    script:
    """
    trimmomatic PE -threads ${task.cpus} -phred33 \\
      -summary ${sample}_summary.txt -compressLevel 5 \\
      ${r1} ${r2} \\
      ${sample}_filtered_1.fastq.gz ${sample}_unpaired_1.fastq.gz \\
      ${sample}_filtered_2.fastq.gz ${sample}_unpaired_2.fastq.gz \\
      ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
