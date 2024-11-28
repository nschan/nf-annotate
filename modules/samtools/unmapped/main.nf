process SAMTOOLS_GET_UNMAPPED {
  tag "${meta}"
  label 'process_low'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'
    : 'biocontainers/samtools:1.15.1--h1170115_0'}"

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*fastq.gz"), emit: unmapped_fq

  script:
  """
    samtools fastq -f 4 ${bam} \
    -1 ${bam}.unmapped_R1.fastq.gz \
    -2 ${bam}.unmapped_R2.fastq.gz
    """
}
