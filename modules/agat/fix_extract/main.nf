process AGAT_FIX_EXTRACT_TRANSCRIPTS {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

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
