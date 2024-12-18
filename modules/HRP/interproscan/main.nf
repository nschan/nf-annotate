process INTERPROSCAN_PFAM {
  tag "$meta"
  label 'process_high'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'docker://interpro/interproscan:5.67-99.0'
    : 'interpro/interproscan:5.67-99.0'}"
  
  publishDir(
    path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
    mode: 'copy',
    overwrite: true,
    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
  ) 
  input:
      tuple val(meta), path(protein_fasta)
  
  output:
      tuple val(meta), path("*.tsv"), emit: protein_tsv
      tuple val(meta), path("*.bed"), emit: nb_bed
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  /opt/interproscan/interproscan.sh \\
     -f TSV,GFF3 \\
     -appl Pfam \\
     -cpu $task.cpus \\
     -i ${protein_fasta} \\
     -b ${prefix}_proteins \\
     -T "${PWD}/tmp"
  grep NB-ARC ${prefix}_proteins.tsv | cut -f1,7,8 > ${prefix}_NB.bed
  """
}

process INTERPROSCAN {
  tag "$meta"
  label 'process_high'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'docker://interpro/interproscan:5.67-99.0'
    : 'interpro/interproscan:5.67-99.0'}"
  input:
      tuple val(meta), path(proteins)
  
  output:
      tuple val(meta), path("*.tsv"), emit: protein_tsv
      tuple val(meta), path("*.gff3"), emit: protein_gff
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  /opt/interproscan/interproscan.sh \\
     -f TSV,GFF3 \\
     -exclappl AntiFam \\
     -cpu $task.cpus \\
     -i ${proteins} \\
     -b ${prefix}_interpro   \\
    -T "${PWD}/tmp"
  """
}

process INTERPROSCAN_EXTENDED {
  tag "$meta"
  label 'process_high'
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'docker://interpro/interproscan:5.67-99.0'
    : 'interpro/interproscan:5.67-99.0'}"
  input:
      tuple val(meta), path(candidate_nb_lrrs)
  
  output:
      tuple val(meta), path("*.tsv"), emit: protein_tsv
      tuple val(meta), path("*.bed"), emit: nb_bed
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  /opt/interproscan/interproscan.sh \\
     -app TIGRFAM,SFLD,SUPERFAMILY,PANTHER,Gene3D,Hamap,ProSiteProfiles,SMART,CDD,PRINTS,PIRSR,Pfam \\
     -f TSV,GFF3 \\
     -cpu $task.cpus \\
     -i ${candidate_nb_lrrs} \\
     -b ${prefix}_proteins   \\
    -T "${PWD}/tmp"
    grep NB-ARC ${prefix}_proteins.tsv | cut -f1,7,8 > ${prefix}_NB.bed
  """
}
process INTERPROSCAN_SUPERFAMILY {
  tag "$meta"
  label 'process_high'
  
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'docker://interpro/interproscan:5.67-99.0'
    : 'interpro/interproscan:5.67-99.0'}"
  input:
      tuple val(meta), path(protein_fasta)
  
  output:
      tuple val(meta), path("*.tsv"), emit: protein_tsv
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  /opt/interproscan/interproscan.sh \\
    -f TSV \\
    -cpu $task.cpus \\
    -app SUPERFAMILY \\
    -i ${protein_fasta} \\
    -b ${prefix}_superfam \\
    -T "${PWD}/tmp"
  """
}