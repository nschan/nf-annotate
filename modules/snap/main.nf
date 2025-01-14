process SNAP {
    tag "$meta"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/snap:2013_11_29--h470a237_1'
    : 'quay.io/biocontainers/snap:2013_11_29--h470a237_1'}"
    input:
        tuple val(meta), path(accession_genome)
        val (organism)

    output:
        tuple val(meta), path("*_snap.gff"), emit: snap_gff
        tuple val(meta), path("*.fna"), emit: snap_fna

    script:
        """
        snap ${organism}.hmm ${accession_genome} \\
            -gff \\
            -aa ${meta}.fna > ${meta}_snap.tmp.gff
        cat ${meta}_snap.tmp.gff \\
            | awk 'BEGIN {OFS="\\t"}; /^[^#]/ {\$2 = "SNAP"; nine=\$9
                    for(i=10; i <= NF; i++) nine=nine" "\$i
                    print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' > ${meta}_snap.gff
        rm ${meta}_snap.tmp.gff
        """
}