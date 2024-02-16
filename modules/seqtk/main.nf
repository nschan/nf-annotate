include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQTK_SUBSET_INPUT {
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
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
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
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
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
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