process AGAT_COMPLEMENT {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(reference_gff), path(additional_gff1, name: "additional_annotations_1.gff3"), path(additional_gff2, name: "additional_annotations_2.gff3"), path(additional_gff3, name: "additional_annotations_3.gff3")

  output:
  tuple val(meta), path("*_merged.gff3"), emit: updated_gff

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  def addfiles = additional_gff2 ? (additional_gff3 ? "--gff ${additional_gff1} --gff ${additional_gff2} --gff ${additional_gff3}" : "--gff ${additional_gff1} --gff ${additional_gff2}") : "-gff ${additional_gff1}"
  """
  agat_sp_merge_annotations.pl \\
    --gff ${reference_gff} \\
    $addfiles \\
    --out ${meta}_merged.gff3 \\
  """
}