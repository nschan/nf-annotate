process MAKEBLASTDB {
    tag "$meta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0' :
        'quay.io/biocontainers/blast:2.14.1--pl5321h6f7f691_0' }"
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta}"), emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi
    sed 's/>\\(AT[A-Z0-9]*\\)/>\\1 unknown GN=\\1 PE=1 SV=1/g' ${fasta_name} \\
      | sed 's/gene=ID=\\(.*\\)/unknown GN=\\1 PE=1 SV=1/g' \\
      | sed 's/>ID=/>/g' > modified_${fasta_name}

    makeblastdb \\
        -in modified_${fasta_name} \\
        -dbtype prot \\
        ${args}

    mkdir -p ${prefix}
    mv modified_${fasta_name}* ${prefix}
    """
}