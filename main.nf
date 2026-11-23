include { DOWNLOAD as DOWNLOAD_READS     } from './modules/download.nf'
include { DOWNLOAD as DOWNLOAD_REFERENCE } from './modules/download.nf'
include { HISAT2_INDEX                   } from './modules/hisat2_index.nf'
include { FASTQC                         } from './modules/fastqc.nf'
include { FASTP                          } from './modules/fastp.nf'
include { HISAT2_ALIGN                   } from './modules/hisat2_align.nf'
include { FEATURECOUNTS                  } from './modules/featurecounts.nf'
include { DESEQ2                         } from './modules/deseq2.nf'
include { MULTIQC                        } from './modules/multiqc.nf'
include { CHECK_RESULTS                  } from './modules/check_results.nf'

params {
    // Table describing the samples: name, condition, checksum and download link
    samplesheet: Path = "${projectDir}/samplesheet.csv"

    // Reference genome and annotation: Ensembl release 115, assembly R64-1-1
    genome_url: String = "https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
    genome_md5: String = "5ae49d2d22cd977abdd6b054b94c544d"
    gtf_url:    String = "https://ftp.ensembl.org/pub/release-115/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz"
    gtf_md5:    String = "b26b4cb4cdf143f605c6e917d1f37f5f"

    // Library strandedness for read counting: 0 = unstranded, 1 = stranded, 2 = reversely stranded
    strandedness: Integer = 0

    // MultiQC report settings
    multiqc_config: Path = "${projectDir}/assets/multiqc_config.yaml"

    // Reference outputs that a finished run is compared against
    expected_dir: Path = "${projectDir}/expected_results"
}

workflow {
    main:
    // --- Setup: reads --------------------------------------------------
    samples  = channel.fromPath(params.samplesheet).splitCsv(header: true)
    reads_in = samples.map { row -> tuple(row.sample, "${row.sample}.fastq.gz", row.fastq_url, row.fastq_md5) }
    reads    = DOWNLOAD_READS(reads_in)

    // --- Setup: reference ----------------------------------------------
    reference_in = channel.of(
        tuple('genome',     'genome.fa.gz',      params.genome_url, params.genome_md5),
        tuple('annotation', 'annotation.gtf.gz', params.gtf_url,    params.gtf_md5)
    )
    reference = DOWNLOAD_REFERENCE(reference_in)
    genome    = reference.filter { id, f -> id == 'genome' }.map { id, f -> f }
    gtf       = reference.filter { id, f -> id == 'annotation' }.map { id, f -> f }

    // --- Setup: aligner index ------------------------------------------
    index = HISAT2_INDEX(genome, gtf)

    // --- Analysis stage 1: quality control of the raw reads ------------
    fastqc = FASTQC(reads)

    // --- Analysis stage 2: adapter and quality trimming ----------------
    fastp = FASTP(reads)

    // --- Analysis stage 3: align to the genome, sort and index ---------
    align = HISAT2_ALIGN(fastp.reads, index.first())

    // --- Analysis stage 4: count reads per gene ------------------------
    bams   = align.bam.map { id, bam, bai -> bam }.collect()
    counts = FEATURECOUNTS(bams, gtf.first(), params.strandedness)

    // --- Analysis stage 5: differential expression with DESeq2 ---------
    deseq2 = DESEQ2(counts.matrix, params.samplesheet)

    // --- Analysis stage 6: summary report of all QC and alignment logs -
    multiqc_inputs = fastqc.zip
        .mix(fastp.json, align.log, align.flagstat, counts.summary)
        .collect()
    multiqc = MULTIQC(multiqc_inputs, params.multiqc_config)

    // --- Final check: compare this run with the expected results -------
    check_inputs = counts.matrix
        .mix(deseq2.results, deseq2.plots.flatten(), multiqc.report, align.log)
        .collect()
    check = CHECK_RESULTS(check_inputs, params.expected_dir)
    check.report.view { f -> f.text }

    publish:
    fastqc_html    = fastqc.html
    fastp_html     = fastp.html
    bam            = align.bam
    align_logs     = align.log.mix(align.flagstat)
    count_matrix   = counts.matrix
    count_summary  = counts.summary
    de_tables      = deseq2.results.mix(deseq2.normalized, deseq2.summary, deseq2.sessioninfo)
    de_figures     = deseq2.plots
    multiqc_report = multiqc.report
    multiqc_data   = multiqc.data
    check_report   = check.report
}

output {
    fastqc_html    { path 'fastqc' }
    fastp_html     { path 'fastp' }
    bam            { path 'alignments' }
    align_logs     { path 'alignment_logs' }
    count_matrix   { path 'counts' }
    count_summary  { path 'counts' }
    de_tables      { path 'deseq2' }
    de_figures     { path 'deseq2/figures' }
    multiqc_report { path 'multiqc' }
    multiqc_data   { path 'multiqc' }
    check_report   { path '.' }
}
