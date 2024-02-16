

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNAP {
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(meta), path(accession_genome)

    output:
        tuple val(meta), path("*_snap.gff"), emit: snap_gff
        tuple val(meta), path("*.fna"), emit: snap_fna

    script:
        """
        snap A.thaliana.hmm ${accession_genome} \\
        -gff \\
        -aa ${meta}.fna > ${meta}_snap.tmp.gff
        cat ${meta}_snap.tmp.gff \\
        | awk 'BEGIN {OFS="\\t"}; /^[^#]/ {\$2 = "SNAP"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' > ${meta}_snap.gff
        rm ${meta}_snap.tmp.gff
        """
}