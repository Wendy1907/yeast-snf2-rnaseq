process MULTIQC {
    conda "${projectDir}/envs/multiqc.yml"

    input:
    path 'inputs/*'
    path config

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_report_data",  emit: data

    script:
    """
    multiqc --force --config ${config} --filename multiqc_report inputs
    """
}
