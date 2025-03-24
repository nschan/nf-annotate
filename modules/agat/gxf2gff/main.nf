process AGAT_GXF2GFF {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(gff_cds), path(gff_augustus), path(gff_snap), path(gff_pasa)

  output:
  tuple val(meta), path("*_predictions.gff"), emit: gff_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  awk 'BEGIN {OFS="\\t"}; {\$2 = "PASA"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' ${gff_pasa} > ${gff_pasa}.fixed

  cat ${gff_augustus} ${gff_snap} ${gff_cds} ${gff_pasa}.fixed > ${meta}_predictions.tmp
  agat_convert_sp_gxf2gxf.pl \\
  -g ${meta}_predictions.tmp \\
  -o ${meta}_predictions.gff
  """
}
