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

process AGAT_GXF2GFF {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

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
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(genome), path(gtf_bambu)

  output:
  tuple val(meta), path(genome), path("*_bambu.gff"), emit: gff_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_convert_sp_gxf2gxf.pl \\
    -g ${gtf_bambu} \\
    -o ${meta}_bambu.gff
  """
}

process AGAT_FUNCTIONAL_ANNOTATION {
  tag "${meta}"
  label 'process_medium'
  publishDir path: { "${params.out}/${task.process}".replace(':', '/').toLowerCase() }, mode: 'copy', overwrite: true, saveAs: { fn -> fn.substring(fn.lastIndexOf('/') + 1) }
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(gff), path(blast_results), path(interpro_results)
  tuple val(meta2), path(blast_reference)

  output:
  tuple val(meta), path("*_functional.gff"), emit: gff_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  sed 's/>\\(AT[A-Z0-9]*\\)/>\\1 unknown GN=\\1 PE=1 SV=1/g' ${blast_reference} \\
    | sed 's/gene=ID=\\(.*\\)/unknown GN=\\1 PE=1 SV=1/g' \\
    | sed 's/>ID=/>/g' > modified_${blast_reference}

  agat_sp_manage_functional_annotation.pl \\
    -f ${gff} \\
    -b ${blast_results} \\
    --db modified_${blast_reference} \\
    -i ${interpro_results} > ${meta}_functional.agat.out
  head -n -1 ${meta}_functional.agat.out > ${meta}_functional.gff
  sed -i '1,/Formating output to GFF3/d' ${meta}_functional.gff
  """
}

process AGAT_GFF2GTF {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(gff)

  output:
  tuple val(meta), path("*.gtf"), emit: gtf_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_convert_sp_gff2gtf.pl \\
  --gff ${gff} \\
  -o ${gff.baseName}.gtf 
  """
}

process AGAT_GFF2BED {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(gff)

  output:
  tuple val(meta), path("*.bed"), emit: bed_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  """
  agat_convert_sp_gff2bed.pl \\
  --gff ${gff} \\
  -o ${gff.baseName}.bed 
  """
}

process AGAT_COMPLEMENT {
  tag "${meta}"
  label 'process_medium'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://depot.galaxyproject.org/singularity/agat:1.4.1--pl5321hdfd78af_0'
    : 'quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0'}"

  input:
  tuple val(meta), path(reference_gff), path(additional_gff1), path(additional_gff2), path(additional_gff3)

  output:
  tuple val(meta), path("*.bed"), emit: bed_file

  script:
  def prefix = task.ext.prefix ?: "${meta}"
  def addfiles = additional_gff2 ? (additional_gff3 ? "--add ${additional_gff1} --add ${additional_gff2} --add ${additional_gff3}" : "--add ${additional_gff1} --add ${additional_gff2}") : "-add ${additional_gff1}"
  """
  agat_sp_complement_annotations.pl \\
    --ref annotation_ref.gff \\
    $addfiles \\
    --out ${meta}_complemented.gff3 \\
  """
}
