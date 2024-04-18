include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIPROT_HRP {
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }

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