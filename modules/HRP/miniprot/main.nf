process MINIPROT_HRP {
    tag "${meta}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/miniprot:0.12--he4a0461_0'
        : 'quay.io/biocontainers/miniprot:0.12--he4a0461_0'}"

    input:
    tuple val(meta), path(accession_genome), path(proteins)

    output:
    tuple val(meta), path("*_miniprot.gff"), emit: miniprot_nlrs

    script:
    """
    miniprot -t${task.cpus} -d ${meta}.mpi ${accession_genome}
    miniprot -Iut${task.cpus} --gff ${meta}.mpi ${proteins} > ${meta}_nlrs_miniprot.gff
    """
}
