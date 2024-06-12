process SNAP {
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