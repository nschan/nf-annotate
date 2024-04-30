process SEQKIT_GET_LENGTH {
  tag "$meta"
  label 'process_low'
  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
    mode: 'copy',
    overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
  conda "bioconda::seqkit=2.4.0"
  input:
      tuple val(meta), path(fasta_file)
  
  output:
      tuple val(meta), path("*_length.txt"), emit: length_estimates
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  seqkit fx2tab --length --name ${fasta_file} > ${meta}_length.txt
  """
}

