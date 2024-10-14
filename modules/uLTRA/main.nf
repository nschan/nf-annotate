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
        tuple val(meta), path(reads), path(reference)

    output:
        tuple val(meta), path("*.bam"), emit: alignment
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