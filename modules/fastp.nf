process FASTP {
    tag "${sample}"
    label 'process_medium'
    conda "${projectDir}/envs/qc.yml"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.trimmed.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    """
    fastp \\
        --in1 ${reads} \\
        --out1 ${sample}.trimmed.fastq.gz \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html \\
        --report_title "${sample}" \\
        --adapter_sequence AGATCGGAAGAGC \\
        --length_required 25 \\
        --thread ${task.cpus}
    """
}
