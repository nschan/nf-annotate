process AGAT_FUNCTIONAL_ANNOTATION {
  tag "${meta}"
  label 'process_medium'
  publishDir path: { "${params.out}/${task.process}".replace(':', '/').toLowerCase() }, mode: 'copy', overwrite: true, saveAs: { fn -> fn.substring(fn.lastIndexOf('/') + 1) }
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(gff), path(blast_results), path(interpro_results)
  tuple val(meta2), path(blast_reference)

  output:
  tuple val(meta), path("*_functional.gff"), emit: gff_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  sed 's/>\\(AT[A-Z0-9]*\\)/>\\1 unknown GN=\\1 PE=1 SV=1/g' ${blast_reference} \\
    | sed 's/gene=ID=\\(.*\\)/unknown GN=\\1 PE=1 SV=1/g' \\
    | sed 's/>ID=/>/g' > modified_${blast_reference}

  agat_sp_manage_functional_annotation.pl \\
    -f ${gff} \\
    -b ${blast_results} \\
    --db modified_${blast_reference} \\
    -i ${interpro_results} > ${meta}_functional.agat.out
  head -n -1 ${meta}_functional.agat.out > ${meta}_functional.gff
  sed -i '1,/Formating output to GFF3/d' ${meta}_functional.gff
  """
}
