process DOWNLOAD_FASTQ {
    tag "${sample}"
    maxForks 3
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), val(url), val(md5)

    output:
    tuple val(sample), path("${sample}.fastq.gz")

    script:
    def ftp_url = url.replaceFirst('^https://', 'ftp://')
    """
    rm -f ${sample}.fastq.gz
    fetch() {
        curl -L --fail --silent --show-error \\
             --retry 5 --retry-all-errors --retry-delay 10 \\
             -o ${sample}.fastq.gz "\$1"
    }
    fetch "${url}" || { echo "Primary URL failed; trying ${ftp_url}" >&2; fetch "${ftp_url}"; }
    echo "${md5}  ${sample}.fastq.gz" | md5sum -c -
    """
}
