process PASA_SEQCLEAN {
    tag "${meta}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://pasapipeline/pasapipeline:2.5.3'
        : 'pasapipeline/pasapipeline:2.5.3'}"

    input:
    tuple val(meta), path(accession_genome), path(accession_transcripts)

    output:
    tuple val(meta), path(accession_genome), path(accession_transcripts), path(".clean")

    script:
    """
    export USER=${workflow.userName}
    /usr/local/src/PASApipeline/bin/seqclean ${accession_transcripts} \\
    -c ${task.cpus}
    """
}

process PASA_TD {
    tag "${meta}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://pasapipeline/pasapipeline:2.5.3'
        : 'pasapipeline/pasapipeline:2.5.3'}"

    input:
    tuple val(meta), path(pasa_assembly_fasta), path(pasa_assembly_gff)

    output:
    tuple val(meta), path("*transdecoder.cds"), emit: transcript_cds
    tuple val(meta), path("*transdecoder.pep"), emit: transcript_pep
    tuple val(meta), path("*transdecoder.gff3"), emit: transcript_gff
    tuple val(meta), path("*transdecoder.bed"), emit: transcript_bed
    tuple val(meta), path("*transdecoder.genome.gff3"), emit: genome_gff
    tuple val(meta), path("*transdecoder.genome.bed"), emit: genome_bed

    script:
    """
    export USER=${workflow.userName}
    /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi \
       --pasa_transcripts_fasta ${pasa_assembly_fasta} \
       --pasa_transcripts_gff3 ${pasa_assembly_gff}
    """
}

process PASA_PIPELINE {
    tag "${meta}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://pasapipeline/pasapipeline:2.5.3'
        : 'pasapipeline/pasapipeline:2.5.3'}"

    input:
    tuple val(meta), path(accession_genome), path(accession_transcripts)

    output:
    tuple val(meta), path("*assemblies.fasta"), emit: pasa_assembly_fasta
    tuple val(meta), path("*pasa_assemblies.gff3"), emit: pasa_assembly_gff
    tuple val(meta), path("*pasa_assemblies.gtf"), emit: pasa_assembly_gtf
    tuple val(meta), path("*pasa_assemblies.bed"), emit: pasa_assembly_bed
    tuple val(meta), path("*ascii_illustrations.out"), emit: pasa_ascii
    tuple val(meta), path("*assemblies_described.txt"), emit: pasa_tsv
    tuple val(meta), path("pasa_DB_*.sqlite"), emit: database

    script:
    def args = task.ext.args ?: ''
    def pasa_config = file("${projectDir}/assets/pasa.config", checkIfExists: true)
    pasa_assemblies_fasta = "pasa_DB_${meta}.sqlite.assemblies.fasta"
    pasa_assemblies_gff = "pasa_DB_${meta}.sqlite.pasa_assemblies.gff3"
    db_name = "pasa_DB_${meta}.sqlite"
    """
        # Clean fasta file, remove empty entries
        cat ${accession_transcripts} \\
         | sed 's/>"/>/g' \\
         | sed '/^>/s/,.*//' \\
         | awk 'BEGIN {RS = ">" ; FS = "\\n" ; ORS = ""} \$2 {print ">"\$0}' > ${accession_transcripts}.tmp
        make_pasa_config.pl --infile ${pasa_config} --trunk ${meta} --outfile pasa_DB.config
        /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \\
           -c pasa_DB.config \\
           -C \\
           -R \\
           -g ${accession_genome} \\
           -u ${accession_transcripts}.tmp \\
           -t ${accession_transcripts}.tmp \\
           --ALIGNERS gmap,minimap2 \\
           --CPU ${task.cpus}
    """
}

process PASA_UPDATE {
    tag "${meta}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://pasapipeline/pasapipeline:2.5.3'
        : 'pasapipeline/pasapipeline:2.5.3'}"
    cache 'lenient'

    input:
    tuple val(meta), path(accession_genome), path(accession_transcripts), path(annotations_gff, name: "input_annotations.gff3"), path(database)

    output:
    tuple val(meta), path("*_gene_structures_updated.gff3"), emit: updated_annotations
    tuple val(meta), path("*_gene_structures_updated.bed"), emit: updated_annotations_bed

    script:
    def pasa_config = file("${projectDir}/assets/pasa.config", checkIfExists: true)
    def db = database
    """
    export USER=${workflow.userName}
    cat /usr/local/src/PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt \\
     | sed 's@<__DATABASE__>@${db}@g' > ${meta}_assembly_conf.config

    make_pasa_config.pl --infile ${pasa_config} --trunk ${meta} --outfile pasa_DB.config

    /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \\
        -c pasa_DB.config \\
        -g ${accession_genome} \\
        -P input_annotations.gff3

    cat /usr/local/src/PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt \\
     | sed 's@<__DATABASE__>@${db}@g' > ${meta}_conf.config
    /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \\
        -c pasa_DB.config \\
        -A \\
        -g ${accession_genome} \\
        -t ${accession_transcripts} \\
        --CPU ${task.cpus}

    cp *.gene_structures_post_PASA_updates.*.gff3 ${meta}_gene_structures_updated.gff3
    cp *.gene_structures_post_PASA_updates.*.bed ${meta}_gene_structures_updated.bed

    """
}