process HISAT2_INDEX {
    tag "genome index"
    label 'process_medium'
    conda "${projectDir}/envs/align.yml"

    input:
    path genome_gz
    path gtf_gz

    output:
    path 'hisat2_index'

    script:
    """
    gunzip -c ${genome_gz} > genome.fa
    gunzip -c ${gtf_gz} > annotation.gtf

    hisat2_extract_splice_sites.py annotation.gtf > splice_sites.txt
    hisat2_extract_exons.py annotation.gtf > exons.txt

    mkdir hisat2_index
    hisat2-build -p ${task.cpus} --ss splice_sites.txt --exon exons.txt \\
                 genome.fa hisat2_index/genome
    """
}
