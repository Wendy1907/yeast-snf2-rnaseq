process FASTQC {
    tag "${sample}"
    conda "${projectDir}/envs/qc.yml"

    input:
    tuple val(sample), path(reads)

    output:
    path "*_fastqc.zip",  emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
