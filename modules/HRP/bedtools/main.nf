include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BEDTOOLS_GETFASTA {
  tag "$meta"
  label 'process_low'
  
  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(protein_fasta), path(bed_file)
  
  output:
      tuple val(meta), path("*.fasta"), emit: fasta_subset
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  bedtools getfasta \\
    -fi ${protein_fasta} \\
    -bed ${bed_file} \\
    -fo ${prefix}_NB_Pfam_Domain_Sequences.fasta
  """
}

process BEDTOOLS_CLUSTER {
  tag "$meta"
  label 'process_low'
  
  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(bed_file)
  
  output:
      tuple val(meta), path("*_clusters"), emit: clusters
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  cat ${bed_file} | sortBed | clusterBed -s | cut -f4,11  > ${bed_file}_clusters
  """
}

process BEDTOOLS_NR_CLUSTERS {
  tag "$meta"
  label 'process_low'
  
  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(clusters), path(length_estimates)
  
  output:
      tuple val(meta), path("*.txt"), emit: r_genes
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """

  join -1 1 -2 1 -o 1.1,1.2,2.5 <( sort -bk1 ${clusters}) <(sort -bk1 ${length_estimates}) | sort -bk2,2 -bk3,3 -nr | sort -uk2,2 | cut -f1 -d ' ' > ${prefix}-R-gene_ID_list.txt
  """
}

