include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

// Dockerfile for minimap2-samtools container is in this folder.

process ALIGN_TO_BAM {
    tag "$meta"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(meta), path(reads), path(reference)

    output:
        tuple val(meta), path("*.bam"), emit: alignment

    script:
        """
        minimap2 -t $task.cpus \\
            -ax map-ont ${reference} ${reads} \\
            | samtools sort -o ${meta}_${reference}.bam
        """
}