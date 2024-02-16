include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IPS2FPG {
  tag "$meta"
  label 'process_low'
  
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/debian:bullseye' :
        'debian:bullseye' }"

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(pfam_out), path(superfamily_out)
  
  output:
      tuple val(meta), path("*out.tsv"), emit: out_tsv
      tuple val(meta), path("*full_length.tsv"), emit: full_length_tsv
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
      def conf1 = file("$projectDir/assets/conf1.tsv", checkIfExists: true)
      def conf2 = file("$projectDir/assets/conf2.tsv", checkIfExists: true)
  """
  cat ${pfam_out} ${superfamily_out} > ${meta}_protein_annotations.tsv
  IPS2fpGs.sh -f -o ${meta}_ips2fp_out.tsv -c ${conf1} -d ${conf2} ${meta}_protein_annotations.tsv
  cat ${meta}_ips2fp_out.tsv | grep full- > ${meta}_ips2fp_out_full_length.tsv
  """
}