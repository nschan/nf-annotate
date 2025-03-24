process AGAT_GFF2BED {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(gff)

  output:
  tuple val(meta), path("*.bed"), emit: bed_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_convert_sp_gff2bed.pl \\
  --gff ${gff} \\
  -o ${gff.baseName}.bed 
  """
}
