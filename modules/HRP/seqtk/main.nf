process SEQTK_SUBSET_RPS {
  tag "$meta"
  label 'process_low'
  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
    mode: 'copy',
    overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
  input:
      tuple val(meta), path(fasta), path(ids1), path(ids2)
  
  output:
      tuple val(meta), path("*_subset.fasta"), emit: fasta_subset

  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
    cat ${ids1} \\
    | cut -f 1 > ${ids1}_gene_ids.txt

    cat ${ids1}_gene_ids.txt ${ids2} > all_ids.txt
    
    seqtk subseq ${fasta} all_ids.txt > ${meta}_subset.fasta
  """
}

process SEQTK_SUBSET_FL {
  tag "$meta"
  label 'process_low'
  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
    mode: 'copy',
    overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
  input:
      tuple val(meta), path(fasta), path(fl_tab)
  
  output:
      tuple val(meta), path("*full_length.fasta"), emit: fasta_subset

  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
    cat ${fl_tab} \\
    | cut -f 1 > ${fl_tab}_gene_ids.txt 

    seqtk subseq ${fasta} ${fl_tab}_gene_ids.txt  > ${meta}_full_length.fasta
  """
}

process SEQTK_SUBSET_CANDIDATES {
  tag "$meta"
  label 'process_low'
  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
    mode: 'copy',
    overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
  input:
      tuple val(meta), path(fasta), path(ids1), path(ids2)
  
  output:
      tuple val(meta), path("*_NBLRR_candidates.fasta"), emit: nblrr_fasta
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
    cat ${ids1} \\
    | cut -f1 -d '-' > ${ids1}_gene_ids.txt

    cat ${ids1}_gene_ids.txt ${ids2} > all_ids.txt

    seqtk subseq ${fasta} all_ids.txt > ${meta}_NBLRR_candidates.fasta
  """
}
