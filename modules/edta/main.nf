process EDTA_FULL {
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
        tuple val(meta), path("*mod.EDTA.intact.gff3") , emit: transposon_intact_gff
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

process LTR {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(genome_fasta)

    output:
        tuple val (meta), path("*.LTR.raw.fa"), path("*.LTR.intact.raw.fa"), path("*.LTR.intact.raw.gff3"), emit: ltrs
    script:
        """
        EDTA_raw.pl \\
            --genome $genome_fasta \\
            --type ltr \\
            --species others \\
            -t $task.cpus
        cp ${genome_fasta}.mod.EDTA.raw/*.fa .
        cp ${genome_fasta}.mod.EDTA.raw/*.gff3 .
        """
}

process TIR {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(genome_fasta)

    output:
        tuple val (meta), path("*.TIR.intact.raw.fa"), path("*.TIR.intact.raw.gff3"), emit: tirs
    script:
        """
        EDTA_raw.pl \\
            --genome $genome_fasta \\
            --type tir \\
            --species others \\
            -t $task.cpus
        cp ${genome_fasta}.mod.EDTA.raw/*.fa .
        cp ${genome_fasta}.mod.EDTA.raw/*.gff3 .
        """
}

process HELITRON {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(genome_fasta)

    output:
        tuple val (meta), path("*.Helitron.intact.raw.fa"), path("*.Helitron.intact.raw.gff3"), emit: helitrons
    script:
        """
        EDTA_raw.pl \\
            --genome $genome_fasta \\
            --type helitron \\
            --species others \\
            -t $task.cpus
        cp ${genome_fasta}.mod.EDTA.raw/*.fa .
        cp ${genome_fasta}.mod.EDTA.raw/*.gff3 .     
        """
}

process LINE {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(genome_fasta)

    output:
        tuple val (meta), path("*LINE.raw.fa"), emit: lines
    script:
        """
        EDTA_raw.pl \\
            --genome $genome_fasta \\
            --type line \\
            --species others \\
            -t $task.cpus
        cp -n ${genome_fasta}.mod.EDTA.raw/LINE/*LINE.raw.fa .
        """
}

process SINE {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(genome_fasta)

    output:
        tuple val (meta), path("*.SINE.raw.fa"), emit: sines
    script:
        """
        EDTA_raw.pl \\
            --genome $genome_fasta \\
            --type sine \\
            --species others \\
            -t $task.cpus
        cp -n ${genome_fasta}.mod.EDTA.raw/SINE/*.SINE.raw.fa . 
        """
}

process MERGE {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        // LTR, TIR, HELITRON, LINE, SINE
        tuple val(meta), path(genome_fasta), path(bed), path(cds_fasta), path(ltr_raw_fa), path(ltr_raw_intact_fa), path(ltr_raw_intact_gff), path(tir_intact_fa), path(tir_intact_gff), path(helitron_intact_fa), path(helitron_intact_gff), path(line_intact_fa), path(sine_intact_fa)

    output:
        tuple val(meta), path("*mod.EDTA.TEanno.gff3"), emit: transposon_annotations
        tuple val(meta), path("*mod.EDTA.TEanno.sum") , emit: transposon_summary
        tuple val(meta), path("*mod.MAKER.masked"), emit: masked
        tuple val(meta), path("*.mod.EDTA.TE.fa.stat.redun.sum"), emit: simple_inconsistency, optional: true
        tuple val(meta), path("*.mod.EDTA.TE.fa.stat.nested.sum"), emit: nested_inconsistency, optional: true
        tuple val(meta), path("*.mod.EDTA.TE.fa.stat.all.sum"), emit: overall_inconsistency, optional: true

    script:
        """
        mkdir -p ${genome_fasta}.mod.EDTA.raw/LINE
        cp ${line_intact_fa} ${genome_fasta}.mod.EDTA.raw/LINE/
        cp ${line_intact_fa} ${genome_fasta}.mod.EDTA.raw/

        mkdir -p ${genome_fasta}.mod.EDTA.raw/SINE
        cp ${sine_intact_fa} ${genome_fasta}.mod.EDTA.raw/SINE/
        cp ${sine_intact_fa} ${genome_fasta}.mod.EDTA.raw/

        mkdir -p ${genome_fasta}.mod.EDTA.raw/LTR
        cp ${ltr_raw_fa} ${ltr_raw_intact_fa} ${ltr_raw_intact_gff} ${genome_fasta}.mod.EDTA.raw/LTR/
        cp ${ltr_raw_fa} ${ltr_raw_intact_fa} ${ltr_raw_intact_gff} ${genome_fasta}.mod.EDTA.raw/

        mkdir -p ${genome_fasta}.mod.EDTA.raw/TIR
        cp ${tir_intact_fa} ${tir_intact_gff} ${genome_fasta}.mod.EDTA.raw/TIR/

        mkdir -p ${genome_fasta}.mod.EDTA.raw/Helitron
        cp ${helitron_intact_fa} ${helitron_intact_gff} ${genome_fasta}.mod.EDTA.raw/Helitron/
        cp ${helitron_intact_fa} ${helitron_intact_gff} ${genome_fasta}.mod.EDTA.raw/

        EDTA.pl \\
            --genome $genome_fasta \\
            --cds $cds_fasta \\
            --exclude $bed \\
            --species others \\
            --anno 1 \\
            --overwrite 0 \\
            --force 1 \\
            -t $task.cpus
        
        mv ${genome_fasta}.mod.EDTA.raw/*EDTA* .
        mv ${genome_fasta}.mod.EDTA.raw/*MAKER* .
        """
}