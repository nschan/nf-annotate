process BAMBU {
    tag "$meta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.4.0--r43hf17093f_0' :
        'quay.io/biocontainers/bioconductor-bambu:3.4.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(fasta), path(gtf), path(bams)

    output:
    tuple val(meta), path("${meta}_bambu_counts_gene.txt")         , emit: ch_gene_counts
    tuple val(meta), path("${meta}_bambu_counts_transcript.txt")   , emit: ch_transcript_counts
    tuple val(meta), path("${meta}_bambu.gtf")                     , emit: extended_gtf
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=$gtf \\
        --fasta=$fasta $bams
    
    mv extended_annotations.gtf ${meta}_bambu.gtf
    mv counts_transcript.txt ${meta}_bambu_counts_transcript.txt
    mv counts_gene.txt ${meta}_bambu_counts_gene.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
        bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
    END_VERSIONS
    """
}