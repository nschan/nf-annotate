process SEQKIT_GET_LENGTH {
  tag "${meta}"
  label 'process_low'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0'
    : 'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0'}"
  conda "bioconda::seqkit=2.4.0"

  input:
  tuple val(meta), path(fasta_file)

  output:
  tuple val(meta), path("*_length.txt"), emit: length_estimates

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  seqkit fx2tab --length --name ${fasta_file} > ${meta}_length.txt
  """
}
