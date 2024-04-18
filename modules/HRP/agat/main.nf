include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AGAT_FILTER_BY_LENGTH {
  tag "$meta"
  label 'process_low'

  conda "bioconda::agat=1.1.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/agat:1.1.0--pl5321hdfd78af_1' :
          'quay.io/biocontainers/agat:1.1.0--pl5321hdfd78af_1' }"

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(gff_file)
  
  output:
      tuple val(meta), path("*_filtered_CDS.gff"), emit: filtered_gff
      tuple val(meta), path("*_filtered_CDS.bed"), emit: filtered_bed
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_sp_filter_gene_by_length.pl \\
  --gff ${gff_file} \\
  --size 12000 --test "<" \\
  -o ${meta}_filtered.gff

  grep CDS ${meta}_filtered.gff > ${meta}_filtered_CDS.gff

  agat_convert_sp_gff2bed.pl \\
  --gff ${meta}_filtered_CDS.gff \\
  -o ${meta}_filtered_CDS.bed
  """
}

process AGAT_EXTRACT_PROTEINS {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(genome_fasta), path(genome_gff)
      val(exclusion_pattern)
  
  output:
      tuple val(meta), path("*_proteins.fasta"), emit: extracted_proteins
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  cat ${genome_gff} | grep CDS | grep -v ${exclusion_pattern} > ${genome_gff}_subset.gff3
  cat ${genome_fasta} | fold > ${genome_fasta.baseName}.fold.fasta
  agat_sp_extract_sequences.pl \\
       -g ${genome_gff}_subset.gff3 \\
       -f ${genome_fasta.baseName}.fold.fasta \\
       -p \\
       --cfs \\
       --cis \\
       -o proteins.fa
  
  grep -B1 "*" proteins.fa | grep -vFf - proteins.fa | fold > ${prefix}_proteins.fasta  
  """
}

process AGAT_EXTRACT_NLR {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(genome_fasta), path(nlr_gff)
  
  output:
      tuple val(meta), path("*NLR_proteins.fasta"), emit: extracted_nlrs
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  cat ${genome_fasta} | fold > ${genome_fasta.baseName}.fold.fasta
  awk -F'\t' -v OFS='\t' '!(\$8=="0,1" || \$8=="0,2") {print \$0}' ${nlr_gff} > filtered_${nlr_gff}
  agat_sp_extract_sequences.pl \\
       -g filtered_${nlr_gff} \\
       -f ${genome_fasta.baseName}.fold.fasta \\
       -p \\
       -t transcript \\
       --cfs \\
       --cis \\
       -o ${prefix}_NLR_proteins.fasta
   """
}

process AGAT_EXTRACT_MINIPROT_NLR {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(genome_fasta), path(nlr_gff)
  
  output:
      tuple val(meta), path("*NLR_proteins.fasta"), emit: extracted_nlrs
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  cat ${genome_fasta} | fold > ${genome_fasta.baseName}.fold.fasta
  agat_sp_extract_sequences.pl \\
       -g ${nlr_gff} \\
       -f ${genome_fasta.baseName}.fold.fasta \\
       -p \\
       -t CDS \\
       --cfs \\
       --cis \\
       -o ${prefix}_NLR_proteins.fasta
   """
}

process AGAT_COMPLEMENT {
  tag "$meta"
  label 'process_low'

  publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
  input:
      tuple val(meta), path(ref_gff), path(nlr_gff)
  
  output:
      tuple val(meta), path("*_liftoff_nlr_merge.gff"), emit: merged_gff
      tuple val(meta), path("*_liftoff_nlr_merge.gtf"), emit: merged_gtf
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_sp_complement_annotations.pl \\
    --ref ${ref_gff} \\
    --add ${nlr_gff} \\
    --out ${meta}_liftoff_nlr_merge.gff  
    agat_convert_sp_gff2gtf.pl \\
    --gff ${meta}_liftoff_nlr_merge.gff   \\
    -o ${meta}_liftoff_nlr_merge.gtf 
  """
}