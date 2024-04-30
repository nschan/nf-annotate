process IPS2FPG {
  tag "$meta"
  label 'process_low'
  
  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
    mode: 'copy',
    overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
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