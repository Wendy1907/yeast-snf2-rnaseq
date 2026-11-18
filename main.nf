include { DOWNLOAD_FASTQ } from './modules/download_fastq.nf'

params {
    // Table describing the samples: name, condition, checksum and download link
    samplesheet: Path = "${projectDir}/samplesheet.csv"
}

workflow {
    main:
    samples = channel.fromPath(params.samplesheet).splitCsv(header: true)
    samples.view { row -> "Sample: ${row.sample} (${row.condition}, replicate ${row.replicate})" }

    to_download = samples.map { row -> tuple(row.sample, row.fastq_url, row.fastq_md5) }
    reads = DOWNLOAD_FASTQ(to_download)
    reads.view { id, fq -> "Downloaded: ${id} -> ${fq.name}" }
}
