process AGAT_GFF2GTF {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(gff)

  output:
  tuple val(meta), path("*.gtf"), emit: gtf_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_convert_sp_gff2gtf.pl \\
  --gff ${gff} \\
  -o ${gff.baseName}.gtf 
  """
}
