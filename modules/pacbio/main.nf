process LIMA {
    tag "$meta"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/lima:2.9.0--h9ee0642_1'
        : 'quay.io/biocontainers/lima:2.9.0--h9ee0642_1'}"
    input:
        tuple val(meta), path(reads)
        path(primers)

    output:
        tuple val(meta), path("*.report") , emit: report
        tuple val(meta), path("*.counts") , emit: counts
        tuple val(meta), path("*.summary"), emit: summary
        tuple val(meta), path("*.bam")    , emit: bam
        tuple val(meta), path("*.bam.pbi")          , optional: true, emit: pbi
        tuple val(meta), path("*.xml")              , optional: true, emit: xml
        tuple val(meta), path("*.json")             , optional: true, emit: json
        tuple val(meta), path("*.clips")            , optional: true, emit: clips
        tuple val(meta), path("*.guess")            , optional: true, emit: guess
    
    script: 
    def args = task.ext.args ?: ''
    if( "$reads" == "${meta}.bam" ) error "Input and output names are the same"

    """
    lima \\ 
        --isoseq \\
        ${reads} \\
        $primers \\
        ${meta}.bam \\
        -j $task.cpus \\
        $args
    """

}