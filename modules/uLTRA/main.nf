process ULTRA_ALIGN {
    tag "$meta"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.1--pyh7cba7a3_1'
        : 'quay.io/biocontainers/ultra_bioinformatics:0.1--pyh7cba7a3_1'}"
    input:
        tuple val(meta), path(reads), path(genome), path(pickle), path(db)
        val(mode)

    output:
        tuple val(meta), path("*.bam"), emit: alignment

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def readmode = mode == "ont" ? '--ont' : '--isoseq'
    """ 
    uLTRA \\
        align \\
        --t $task.cpus \\
        --prefix ${meta} \\
        --index ./ \\
        $args \\
        ${genome} \\
        $readmode
        ${reads} \\
        ./

    samtools \\
        sort \\
        --threads $task.cpus \\
        -o ${meta}.bam \\
        -O BAM \\
        $args2 \\
        ${meta}.sam

    rm ${prefix}.sam
    """
}

process ULTRA_INDEX {
    tag "$meta"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.1--pyh7cba7a3_1'
        : 'quay.io/biocontainers/ultra_bioinformatics:0.1--pyh7cba7a3_1'}"
    input:
        tuple val(meta), path(genome), path(annotations)

    output:
        tuple path("*.pickle"), path("*.db"), emit: index
    script:
    def args = task.ext.args ?: ''
    """
    ultra index \\
        $genome \\
        $annotations \\
        $args
    """
}