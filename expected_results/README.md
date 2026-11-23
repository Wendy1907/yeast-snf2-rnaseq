# Expected results

Reference outputs from a complete run of the default pipeline (`nextflow run main.nf -profile conda`)
on a clean clone of this repository (WSL2 Ubuntu, amd64, Nextflow 26.10.6). The last pipeline step
(`CHECK_RESULTS`) compares a new run against these files and writes `results/check_report.txt`.

| Pipeline stage | Files here | How they are checked |
|---|---|---|
| 3. HISAT2 alignment + samtools | `alignment_logs/*.hisat2.log`, `alignment_logs/*.flagstat.txt` | overall alignment rate within 0.05 percentage points |
| 4. featureCounts | `counts/gene_counts.tsv`, `counts/counts.txt.summary` | same genes and samples; total difference ≤ 0.1% of reads |
| 5. DESeq2 | `deseq2/deseq2_results.csv`, `deseq2/normalized_counts.csv`, `deseq2/deseq2_summary.txt` | results table: fold changes within 0.1%, p-values within 0.01 on the log10 scale |
| 5. DESeq2 figures | `deseq2/figures/*.png` | must exist; compare by eye (PCA, volcano, MA plot, heatmaps) |

Not stored here: raw reads and BAM files (large), and the MultiQC report (it contains run dates and paths;
the check only requires that it is produced).

Why a tolerance: integer read counts should reproduce exactly, but the statistical output (DESeq2) can differ in
the last decimal places between CPU architectures, for example Intel/AMD and Apple silicon.
