process SEQTK_SUBSET_INPUT {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1'
    : 'quay.io/biocontainers/seqtk:1.4--he4a0461_1'}"

  input:
  tuple val(meta), path(genome_fasta), path(annotation_liftoff), path(annotation_bambu), path(contig_list)

  output:
  tuple val(meta), path("*_subset.fa"), path("*_liftoff_subset.gff"), path("*_bambu_subset.gff"), emit: subset

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
    seqtk subseq ${genome_fasta} ${contig_list} > ${meta}_genome_subset.fa
    grep -wf ${contig_list} ${annotation_liftoff} > ${meta}_liftoff_subset.gff.tmp
    grep -wf ${contig_list} ${annotation_bambu} > ${meta}_bambu_subset.gff
    awk 'BEGIN {OFS="\\t"}; /^[^#]/ {\$2 = "LIFTOFF"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' ${meta}_liftoff_subset.gff.tmp > ${meta}_liftoff_subset.gff
    rm ${meta}_liftoff_subset.gff.tmp
  """
}

process SEQTK_SUBSET_FASTA {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1'
    : 'quay.io/biocontainers/seqtk:1.4--he4a0461_1'}"

  input:
  tuple val(meta), path(genome_fasta), path(contig_list)

  output:
  tuple val(meta), path("*_subset.fa"), emit: subset

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
    seqtk subseq ${genome_fasta} ${contig_list} > ${meta}_genome_subset.fa
  """
}

process SUBSET_ANNOTATIONS {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1'
    : 'quay.io/biocontainers/seqtk:1.4--he4a0461_1'}"

  input:
  tuple val(meta), path(annotation_liftoff), path(contig_list)

  output:
  tuple val(meta), path("*_liftoff_subset.gff"), emit: annotations

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
    grep -wf ${contig_list} ${annotation_liftoff} > ${meta}_liftoff_subset.gff.tmp
    awk 'BEGIN {OFS="\\t"}; /^[^#]/ {\$2 = "LIFTOFF"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' ${meta}_liftoff_subset.gff.tmp > ${meta}_liftoff_subset.gff
    rm ${meta}_liftoff_subset.gff.tmp
  """
}
