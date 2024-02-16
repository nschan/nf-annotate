include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQKIT_GET_LENGTH {
  tag "$meta"
  label 'process_low'

  conda "bioconda::seqkit=2.4.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0' :
          'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0' }"

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
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

