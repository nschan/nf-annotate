include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AGAT_FIX_EXTRACT_TRANSCRIPTS {
  tag "$meta"
  label 'process_medium'

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename,
                                    options:params.options, 
                                    publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                    publish_id:meta) }
  input:
      tuple val(meta), path(genome_fasta), path(genome_gff)
  
  output:
      tuple val(meta), path("*_transcripts.fasta"), emit: extracted_transcripts
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  cat ${genome_fasta} | fold > ${genome_fasta.baseName}.fold.fasta

  agat_convert_sp_gxf2gxf.pl \\
    -g ${genome_gff} \\
    -o ${genome_gff}.tmp

  agat_sp_extract_sequences.pl \\
  -g ${genome_gff}.tmp \\
  -f ${genome_fasta.baseName}.fold.fasta \\
  -o ${meta}_transcripts.fasta \\
  -t exon \\
  --merge
  """
}

process AGAT_GXF2GFF {
  tag "$meta"
  label 'process_medium'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(gff_cds), path(gff_augustus), path(gff_snap), path(gff_pasa)
  
  output:
      tuple val(meta), path("*_predictions.gff"), emit: gff_file
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  awk 'BEGIN {OFS="\\t"}; {\$2 = "PASA"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' ${gff_pasa} > ${gff_pasa}.fixed

  cat ${gff_augustus} ${gff_snap} ${gff_cds} ${gff_pasa}.fixed > ${meta}_predictions.tmp
  agat_convert_sp_gxf2gxf.pl \\
  -g ${meta}_predictions.tmp \\
  -o ${meta}_predictions.gff
  """
}

process AGAT_GTF2GFF {
  tag "$meta"
  label 'process_medium'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(genome), path(gtf_bambu)
  
  output:
      tuple val(meta), path(genome) , path("*_bambu.gff"), emit: gff_file
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_convert_sp_gxf2gxf.pl \\
    -g ${gtf_bambu} \\
    -o ${meta}_bambu.gff
  """
}

process AGAT_FUNCTIONAL_ANNOTATION {
  tag "$meta"
  label 'process_medium'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(gff), path(blast_results), path(interpro_results)
      tuple val(meta2), path(blast_reference)
  
  output:
      tuple val(meta), path("*_functional.gff"), emit: gff_file
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  sed 's/gene=ID=\\(.*\\)/GN=\\1 PE=1 SV=1/g' ${blast_reference} > modified_${blast_reference}
  agat_sp_manage_functional_annotation.pl \\
    -f ${gff} \\
    -b ${blast_results} \\
    --db modified_${blast_reference} \\
    -i ${interpro_results} > ${meta}_functional.agat.out
  sed '1,/Formating output to GFF3/d' ${meta}_functional.agat.out > ${meta}_functional.gff
  """
}