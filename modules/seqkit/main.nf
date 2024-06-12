process SEQKIT_GET_LENGTH {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
      tuple val(meta), path(genome_fasta)
      val(length)
  
    output:
      tuple val(meta), path(genome_fasta), path("*_subset.txt"), emit: large_contigs
      tuple val(meta), path("*_subset.txt"), emit: contig_list
    script:
      def prefix = task.ext.prefix ?: "${meta}"
      
  """
  seqkit fx2tab --length --name ${genome_fasta} | awk '\$2 > ${length} {print \$1 "\\t " \$2}' | cut -f1 > ${meta}_subset.txt
  """
}