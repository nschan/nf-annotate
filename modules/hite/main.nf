process HITE {
  tag "${meta}"
  label 'process_high'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'docker://kanghu/hite:3.2.0'
    : 'kanghu/hite:3.2.0'}"

  input:
  tuple val(meta), path(genome_fasta)

  output:
  tuple val(meta), path("longest_repeats_*.fa"), emit: longest_repeats
  tuple val(meta), path("confident_tir_*.fa"), emit: confident_tir
  tuple val(meta), path("confident_helitron_*.fa"), emit: confident_helitron
  tuple val(meta), path("confident_non_ltr_*.fa"), emit: confident_non_ltr
  tuple val(meta), path("confident_other_*.fa"), emit: confident_other
  tuple val(meta), path("confident_ltr_cut.fa.cons"), emit: confident_ltr_cut_cons
  tuple val(meta), path("*HiTE.out"), emit: hite_out
  tuple val(meta), path("*HiTE.gff"), emit: hite_gff
  tuple val(meta), path("*HiTE.tbl"), emit: hite_tbl

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  export PATH=/HiTE/tools:/HiTE/module:/opt/conda/envs/HiTE/bin:\$PATH
  python /HiTE/main.py \\
    --genome ${genome_fasta} \\
    --thread ${task.cpus} \\
    --outdir \$PWD/${prefix} \\
    --annotate 1
  wait
  mv ${prefix}/longest_repeats_*.fa .
  if ls ${prefix}/confident_tir_*.fa 1> /dev/null 2>&1; then
    mv ${prefix}/confident_tir_*.fa .
  fi
  if ls ${prefix}/confident_helitron_*.fa 1> /dev/null 2>&1; then
    mv ${prefix}/confident_helitron_*.fa .
  fi
  if ls ${prefix}/confident_non_ltr_*.fa 1> /dev/null 2>&1; then
    mv ${prefix}/confident_non_ltr_*.fa .
  fi
  if ls ${prefix}/confident_ltr_cut.fa.cons 1> /dev/null 2>&1; then
    mv ${prefix}/confident_ltr_cut.fa.cons .
  fi
  if ls ${prefix}/confident_other_*.fa 1> /dev/null 2>&1; then
    mv ${prefix}/confident_other_*.fa .
  fi
  mv ${prefix}/HiTE.out ${prefix}.HiTE.out
  mv ${prefix}/HiTE.gff ${prefix}.HiTE.gff
  mv ${prefix}/HiTE.tbl ${prefix}.HiTE.tbl
  """
}