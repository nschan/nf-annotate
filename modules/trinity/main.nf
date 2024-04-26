include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

// Dockerfile for minimap2-samtools container is in this folder.

process TRINITY {
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.15.1--pl5321h146fbdb_3':
        'quay.io/biocontainers/trinity:2.15.1--pl5321h146fbdb_3' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fa.gz")    , emit: transcript_fasta
    tuple val(meta), path("*.log")      , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"


    // Define the memory requirements. Trinity needs this as an option.
    def avail_mem = 7
    if (!task.memory) {
        log.info '[Trinity] Available memory not known - defaulting to 7GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

    """
    # Note that Trinity needs the word 'trinity' in the outdir

    Trinity \\
        --genome_guided_bam $bam \\
        --genome_guided_max_intron 10000 \\
        --max_memory $avail_mem
        --output ${prefix}_trinity \\
        --CPU $task.cpus \\
        $args \\
        > >(tee ${prefix}.log)
    
    gzip \\
        -cf \\
        ${prefix}_trinity.Trinity.fasta \\
        > ${prefix}.fa.gz

    rm ${prefix}_trinity.Trinity.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: \$(Trinity --version | grep 'Trinity version:' | sed 's/Trinity version: Trinity-//')
    END_VERSIONS
    """
}