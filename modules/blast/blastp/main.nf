process BLASTP {
    tag "$meta"
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    val out_ext

    output:
    tuple val(meta), path("*.xml"), optional: true, emit: xml
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    tuple val(meta), path("*.csv"), optional: true, emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    switch ( out_ext ) {
        case "xml": outfmt = 5; break
        case "tsv": outfmt = 6; break
        case "csv": outfmt = 10; break
        default:
            outfmt = '6';
            out_ext = 'tsv';
            log.warn("Unknown output file format provided (${out_ext}): selecting BLAST default of tabular BLAST output (tsv)");
            break
    }

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi
    db_path=${db}
    DB=`find -L \$db_path -name "*.phr" | sed 's/\\.phr\$//'`
    
    blastp \\
        -query ${fasta_name} \\
        -out ${prefix}_blast.${out_ext} \\
        -db \$DB \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 10 \\
        -outfmt ${outfmt} \\
        ${args}
    """
}