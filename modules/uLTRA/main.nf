process ULTRA_ALIGN {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
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
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
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