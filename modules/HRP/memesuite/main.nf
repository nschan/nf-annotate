include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEME {
  tag "$meta"
  label 'process_medium'
  
  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(protein_fasta)
  
  output:
      tuple val(meta), path("*.txt"), emit: meme_out
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  meme ${protein_fasta}  \\
     -protein  \\
     -o ${prefix} \\
     -mod zoops \\
     -nmotifs 19 \\
     -minw 4 \\
     -maxw 7 \\
     -objfun classic \\
     -markov_order 0 \\
     -p $task.cpus
  
  cp ${prefix}/meme.txt ${prefix}_meme_out.txt 
  """
}

process MAST {
  tag "$meta"
  label 'process_medium'
  
  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(protein_fasta), path(meme_out)
  
  output:
      tuple val(meta), path("*mast_out.txt"), emit: mast_out
      tuple val(meta), path("*mast_geneIDs.txt"), emit: mast_geneids
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  mast -o ${meta}_mast ${meme_out} ${protein_fasta}
  cp ${meta}_mast/mast.txt ${meta}_mast_out.txt
  cat ${meta}_mast_out.txt | grep -oE "AT[1-5C]G[0-9]+.[0-9]+|evm[0-9a-z\\.]*" > ${meta}_mast_geneIDs.txt
  """
}