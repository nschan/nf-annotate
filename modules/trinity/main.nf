process TRINITY {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_transcripts.fa.gz") , emit: transcript_fasta
    tuple val(meta), path("*.log")      , emit: log

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
        --max_memory ${avail_mem}G \\
        --output ${prefix}_trinity \\
        --CPU $task.cpus \\
        $args \\
        > >(tee ${prefix}.log)
    
    cp "${prefix}_trinity/Trinity-GG.fasta" ${prefix}_transcripts.fa
    """
}