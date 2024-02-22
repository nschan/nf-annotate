include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENBLAST_G {
  tag "$meta"
  label 'process_medium'
  container "gitlab.lrz.de:5005/beckerlab/container-playground/genblastg:2a25a73b"
  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(genome_fasta), path(nb_lrr_fasta)
  
  output:
      tuple val(meta), path("*genblastG-output*.pro"), emit: genblast_pro
      tuple val(meta), path("*genblastG-output*.gff"), emit: genblast_gff
      
  /* genblastG
   Dockerfile is included in the module directory.
     !! This script is not very elegant !!
     The blast DB generated by genblastG (via formatdb) filters
     letters valid for proteins with no obvious pattern 
     (E,F,Y,Q seem to cause problems).
     Therefore, blast+ makeblastdb is used to generate the DB.
   formatdb needs to be in the current directory.
     Linking genblast_extension to . is not optimal,
     but genblast needs formatdb in the current directory
     since it calls ./formatdb;
     adding the genblastG folder to PATH does not help
  */

  script:
      def prefix = task.ext.prefix ?: "${meta}"

  """
  export PATH="/opt/bin/ncbi-blast/:$PATH"
  ln -s /opt/genblastG_extension/* .
  makeblastdb -dbtype prot -in ${nb_lrr_fasta} -out ${nb_lrr_fasta}
  ./genblastG -q ${nb_lrr_fasta} -t ${genome_fasta} -gff -cdna -pro -o ${prefix}_genblastG-output
  """
}