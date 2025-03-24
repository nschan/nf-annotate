process AGAT_GTF2GFF {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(genome), path(gtf_bambu)

  output:
  tuple val(meta), path(genome), path("*_bambu.gff"), emit: gff_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_convert_sp_gxf2gxf.pl \\
    -g ${gtf_bambu} \\
    -o ${meta}_bambu.gff
  """
}