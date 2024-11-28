process SAMTOOLS_FASTQ {
    tag "$meta"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0':
        'biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.fq.gz"), emit: fasta

    script:
    """
    samtools fastq $bam | gzip | ${meta}.fq.gz
    """
}
