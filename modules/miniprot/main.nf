include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIPROT {
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
        tuple val(meta), path("*.gff"), emit: miniprot_gff

    script:
        """
        miniprot -t $task.cpus --gff-only ${accession_genome} $params.reference_proteins > ${meta}_miniprot.tmp.gff
        awk 'BEGIN {OFS="\\t"}; /^[^#]/ {\$2 = "MINIPROT"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' ${meta}_miniprot.tmp.gff > ${meta}_miniprot.gff
        rm ${meta}_miniprot.tmp.gff
        """
}