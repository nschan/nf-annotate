process REFINE {
    tag "$meta"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/isoseq3:4.0.0--h9ee0642_0'
        : 'quay.io/biocontainers/isoseq3:4.0.0--h9ee0642_0'}"
    input:
        tuple val(meta), path(bamfile)
    output:
        tuple val(meta), path("*.flnc.bam") , emit: refined_bam
    script:
    def prefix = task.ext.prefix ?: "${meta}"
    def min_polya = params.refine_polya ?: 20
    def primers = params.primers

    """
    isoseq refine \\
        -j $task.cpus \\
        ${bamfile} \\
        ${primers} \\
        ${prefix}.flnc.bam \\
        --require-polya \\
        --min-polya-length ${min_polya}
    """
}