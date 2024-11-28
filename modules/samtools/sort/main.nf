process SAMTOOLS_SORT {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0'
    : 'biocontainers/samtools:1.15.1--h1170115_0'}"

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*.bam"), emit: bam

  script:
  def prefix = task.ext.prefix ? "${meta}${task.ext.prefix}" : "${meta}"
  """
    samtools sort ${options.args} -@ ${task.cpus} -o ${prefix}.bam -T ${prefix} ${bam}
    """
}
