process TRASH {
    tag "${meta}"
    label 'process_high'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? ''
        : 'quay.io/schandry_containers/trash:latest'}"

    input:
    tuple val(meta), path(genome_fasta)

    output:
    tuple val(meta), path("*.all_repeats.fa.csv"), emit: all_repeats_fa, optional: true
    tuple val(meta), path("*_TRASH.gff"),          emit: repeats_gff,    optional: true
    tuple val(meta), path("*_TRASH_circos.pdf"),   emit: circos_plot,    optional: true
    tuple val(meta), path("*_TRASH_summary.csv"),  emit: summary,        optional: true


    script:
    def prefix = task.ext.prefix ?: "${meta}"
    """
    mkdir ${prefix}_tmp
    TRASH_run.sh ${genome_fasta} \\
        --par ${task.cpus} \\
        --o \$NXF_TASK_WORKDIR/${prefix}_tmp

    mv ${prefix}_tmp/all.repeats.from*.csv ${prefix}_TRASH.all_repeats.fa.csv
    mv ${prefix}_tmp/TRASH_*.gff ${prefix}_TRASH.gff
    mv ${prefix}_tmp/plots/*_circos.pdf ${prefix}_TRASH_circos.pdf
    mv ${prefix}_tmp/Summary.of.repetitive.regions_*.csv ${prefix}_TRASH_summary.csv
  """
}
