process MINIPROT {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 

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