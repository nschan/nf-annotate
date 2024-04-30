process PASA_SEQCLEAN {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:  
        tuple val(meta), path(accession_genome), path(accession_transcripts)

    output:
        tuple val(meta), path(accession_genome), path(accession_transcripts), path(".clean")
/*
Set $USER https://github.com/PASApipeline/PASApipeline/issues/155
*/
    script: 
    """
    export USER=${workflow.userName}
    /usr/local/src/PASApipeline/bin/seqclean ${accession_transcripts} \\
    -c $task.cpus
    """
}

process PASA_TD {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:  
        tuple val(meta), path(pasa_assembly_fasta), path(pasa_assembly_gff)

    output:
        tuple val(meta), path("*transdecoder.cds"), emit: transcript_cds
        tuple val(meta), path("*transdecoder.pep"), emit: transcript_pep
        tuple val(meta), path("*transdecoder.gff3"), emit: transcript_gff
        tuple val(meta), path("*transdecoder.bed"), emit: transcript_bed
        tuple val(meta), path("*transdecoder.genome.gff3"), emit: genome_gff
        tuple val(meta), path("*transdecoder.genome.bed"), emit: genome_bed
/*
Set $USER https://github.com/PASApipeline/PASApipeline/issues/155
*/
    script: 
    """
    export USER=${workflow.userName}
    /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi \
       --pasa_transcripts_fasta ${pasa_assembly_fasta} \
       --pasa_transcripts_gff3 ${pasa_assembly_gff}
    """
}

process PASA_PIPELINE {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 

    input:
        tuple val(meta), path(accession_genome), path(accession_transcripts)

    output:
        tuple val(meta), path("*assemblies.fasta"), emit: pasa_assembly_fasta
        tuple val(meta), path("*pasa_assemblies.gff3"), emit: pasa_assembly_gff
        tuple val(meta), path("*pasa_assemblies.gtf"), emit: pasa_assembly_gtf
        tuple val(meta), path("*pasa_assemblies.bed"), emit: pasa_assembly_bed
        tuple val(meta), path("*ascii_illustrations.out"), emit: pasa_ascii
        tuple val(meta), path("*assemblies_described.txt"), emit: pasa_tsv

    script:
    def args = task.ext.args ?: ''
    def pasa_config = file("$projectDir/assets/pasa.config", checkIfExists: true)
    pasa_assemblies_fasta = "pasa_DB_${meta}.sqlite.assemblies.fasta"
    pasa_assemblies_gff = "pasa_DB_${meta}.sqlite.pasa_assemblies.gff3"
    db_name = "pasa_DB_${meta}.sqlite"
        """
        # Clean fasta file, remove empty entries
        cat ${accession_transcripts} \\
         | sed 's/>"/>/g' \\
         | sed '/^>/s/,.*//' \\
         | awk 'BEGIN {RS = ">" ; FS = "\\n" ; ORS = ""} \$2 {print ">"\$0}' > ${accession_transcripts}.tmp
        make_pasa_config.pl --infile $pasa_config --trunk ${meta} --outfile pasa_DB.config
        /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \\
           -c pasa_DB.config \\
           -C \\
           -R \\
           -g ${accession_genome} \\
           -u ${accession_transcripts}.tmp \\
           -t ${accession_transcripts}.tmp \\
           --ALIGNERS gmap,minimap2 \\
           --CPU $task.cpus
        """
}