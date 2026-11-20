process FEATURECOUNTS {
    label 'process_medium'
    conda "${projectDir}/envs/count.yml"

    input:
    path bams
    path gtf_gz
    val strandedness

    output:
    path "counts.txt",         emit: counts
    path "counts.txt.summary", emit: summary
    path "gene_counts.tsv",    emit: matrix

    script:
    """
    gunzip -c ${gtf_gz} > annotation.gtf

    featureCounts \\
        -T ${task.cpus} \\
        -s ${strandedness} \\
        -t exon -g gene_id \\
        -a annotation.gtf \\
        -o counts.txt \\
        *.bam

    grep -v '^#' counts.txt | cut -f1,7- | sed '1s/\\.sorted\\.bam//g' > gene_counts.tsv
    """
}
