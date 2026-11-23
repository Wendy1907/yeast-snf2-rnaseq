process CHECK_RESULTS {
    conda "${projectDir}/envs/check.yml"

    input:
    path 'observed/*'
    path expected

    output:
    path "check_report.txt", emit: report

    script:
    """
    check_results.py --expected ${expected} --observed observed --out check_report.txt
    """
}
