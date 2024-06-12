// Dockerfile for minimap2-samtools container is in this folder.

process ALIGN_TO_BAM {
    tag "$meta"
    label 'process_low'
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

    script:
        """
        minimap2 -t $task.cpus \\
            -ax splice:hq -uf ${reference} ${reads} \\
            | samtools sort -o ${meta}_${reference}.bam
        """
}