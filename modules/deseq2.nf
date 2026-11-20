process DESEQ2 {
    label 'process_medium'
    conda "${projectDir}/envs/deseq2.yml"

    input:
    path counts
    path samplesheet

    output:
    path "deseq2_results.csv",     emit: results
    path "normalized_counts.csv",  emit: normalized
    path "*.png",                  emit: plots
    path "deseq2_summary.txt",     emit: summary
    path "deseq2_sessioninfo.txt", emit: sessioninfo

    script:
    """
    deseq2.R ${counts} ${samplesheet}
    """
}
