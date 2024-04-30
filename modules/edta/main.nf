process EDTA {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(genome_fasta), path(bed), path(cds_fasta)

    output:
        tuple val(meta), path("*mod.EDTA.TEanno.gff3"), emit: transposon_annotations
        tuple val(meta), path("*mod.EDTA.TEanno.sum") , emit: transposon_summary
        tuple val(meta), path("*mod.MAKER.masked"), emit: masked
        tuple val(meta), path("*.mod.EDTA.TE.fa.stat.redun.sum"), emit: simple_inconsistency, optional: true
        tuple val(meta), path("*.mod.EDTA.TE.fa.stat.nested.sum"), emit: nested_inconsistency, optional: true
        tuple val(meta), path("*.mod.EDTA.TE.fa.stat.all.sum"), emit: overall_inconsistency, optional: true

    script:
        """
        EDTA.pl \\
            --genome $genome_fasta \\
            --cds $cds_fasta \\
            --exclude $bed \\
            --species others \\
            --anno 1 \\
            --force 1 \\
            -t $task.cpus
        mv ${genome_fasta}.mod.EDTA.raw/* .
        """
}