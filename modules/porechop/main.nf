process PORECHOP {
    tag "${meta}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/schandry_containers/porechop-pigz:latest'
        : 'quay.io/schandry_containers/porechop-pigz:latest'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_porechopped.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    ## To ensure ID matches rest of pipeline based on meta.id rather than input file name
    
    [[ -f ${prefix}.fastq.gz   ]] || ln -s ${reads} ${prefix}.fastq.gz

    porechop \\
        -i ${prefix}.fastq.gz \\
        -t ${task.cpus} \\
        --no_split \\
        ${args} \\
        -o ${prefix}_porechopped.fastq.gz \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$( porechop --version )
    END_VERSIONS
    """
}
