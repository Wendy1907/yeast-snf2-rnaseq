process HISAT2_ALIGN {
    tag "${sample}"
    label 'process_medium'
    conda "${projectDir}/envs/align.yml"

    input:
    tuple val(sample), path(reads)
    path index

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
    path "${sample}.hisat2.log",   emit: log
    path "${sample}.flagstat.txt", emit: flagstat

    script:
    """
    set -o pipefail

    hisat2 -p ${task.cpus} -x ${index}/genome -U ${reads} 2> ${sample}.hisat2.log \\
        | samtools sort -@ ${task.cpus} -m 256M -o ${sample}.sorted.bam -

    samtools index ${sample}.sorted.bam
    samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat.txt
    """
}
