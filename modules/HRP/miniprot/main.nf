process MINIPROT_HRP {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 

    input:
        tuple val(meta), path(accession_genome), path(proteins)

    output:
        tuple val(meta), path("*_miniprot.gff"), emit: miniprot_nlrs

    script:
        """
        miniprot -t$task.cpus -d ${meta}.mpi $accession_genome
        miniprot -Iut$task.cpus --gff ${meta}.mpi $proteins > ${meta}_nlrs_miniprot.gff
        """
}