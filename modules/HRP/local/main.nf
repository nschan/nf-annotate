include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_R_GENE_GFF {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(gff), path(r_gene_list)
  
  output:
      tuple val(meta), path("*R_gene_annotations.gff"), emit: r_gene_gff
      tuple val(meta), path("*_R_genes_merged.gff"),    emit: r_genes_merged_gff
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  /*
  The second part of this script does:
   - Sort gff file for bedtools
   - Merge with bedtools (stranded) and preserve columns, take mean of blast score and collapse annotation column
   - Reorder bed to gff via awk and write to _merged file
  */
  """
  grep -Fw -f ${r_gene_list} ${gff} > ${meta}_R_gene_annotations.gff
  cat ${meta}_R_gene_annotations.gff \\
  | sort -k1,1 -k4,4n \\
  | bedtools merge -s -c 2,3,6,7,8,9, -o distinct,distinct,mean,distinct,distinct,collapse \\
  | awk 'BEGIN {FS="\\t"; OFS="\\t"} { print \$1,\$4,\$5,\$2,\$3,\$6,\$7,\$8,\$9}' > ${meta}_R_genes_merged.gff
  """
}

process FILTER_R_GENES {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(pfam_out), path(superfamily_out)
      def conf1 = file("$projectDir/assets/conf1.tsv", checkIfExists: true)
      def conf2 = file("$projectDir/assets/conf2.tsv", checkIfExists: true)  
  output:
      tuple val(meta), path("*NLR_table.tsv"), emit: out_tsv
      tuple val(meta), path("*NLR_genes.tsv"), emit: full_length_tsv
  script:
      def prefix = task.ext.prefix ?: "${meta}"

  """
  filter_R_genes.R ${pfam_out} ${superfamily_out} ${meta}
  """
}