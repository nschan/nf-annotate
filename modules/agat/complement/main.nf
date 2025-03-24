process AGAT_COMPLEMENT {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(reference_gff), path(additional_gff1), path(additional_gff2), path(additional_gff3)

  output:
  tuple val(meta), path("*.bed"), emit: bed_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  def addfiles = additional_gff2 ? (additional_gff3 ? "--add ${additional_gff1} --add ${additional_gff2} --add ${additional_gff3}" : "--add ${additional_gff1} --add ${additional_gff2}") : "-add ${additional_gff1}"
  """
  agat_sp_complement_annotations.pl \\
    --ref annotation_ref.gff \\
    $addfiles \\
    --out ${meta}_complemented.gff3 \\
  """
}