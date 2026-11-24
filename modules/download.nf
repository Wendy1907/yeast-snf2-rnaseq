process DOWNLOAD {
    tag "${id}"
    conda "${projectDir}/envs/download.yml"
    maxForks 3
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(id), val(filename), val(url), val(md5)

    output:
    tuple val(id), path(filename)

    script:
    def ftp_url = url.replaceFirst('^https://', 'ftp://')
    """
    rm -f ${filename}
    fetch() {
        curl -L --fail --silent --show-error \\
             --retry 5 --retry-all-errors --retry-delay 10 --speed-limit 100000 --speed-time 30 \\
             -o ${filename} "\$1"
    }
    fetch "${url}" || { echo "Primary URL failed; trying ${ftp_url}" >&2; fetch "${ftp_url}"; }
    echo "${md5}  ${filename}" | md5sum -c -
    """
}
