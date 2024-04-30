/*
This is kind of parallelized via xargs.
*/

process EVIDENCEMODELER_PART_EXEC_MERGE {
    tag "$meta"
    label 'process_high'
    
    conda (params.enable_conda ? "bioconda::evidencemodeler=2.1.0" : null)
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 

    input:
    tuple val(meta), path(genome), path(genes), path(proteins), path(transcripts)

    output:
    tuple val(meta), path("*full.gff"), emit: gff
    tuple val(meta), path("*full.pep"), emit: peptides
    tuple val(meta), path("*full.cds"), emit: cds
    tuple val(meta), path("*full.bed"), emit: bed
    tuple val(meta), path("*_evm.out"), emit: evm_out

    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    partitions = "${meta}_partitions_list.out"
    evm_commands = "${meta}_commands.evm.list"
    evm_out = "${meta}_evm.out"
    protein_options = ""
    transcript_options = ""
    def weights = file("$projectDir/assets/weights.tsv", checkIfExists: true)
    if (proteins) {
       protein_options = "--protein_alignments $proteins"   
    }
    if (transcripts) {
       transcript_options = "--transcript_alignments $transcripts"
    }
    """
    # Convert miniprot gff to evm compatible gff..
    # The container requires python3, not python. Replace shebang and write to local file..
    sed 's/python\$/python3/' /usr/local/bin/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py > miniprot_GFF_2_EVM_GFF3.py
    chmod +x miniprot_GFF_2_EVM_GFF3.py
    ./miniprot_GFF_2_EVM_GFF3.py $proteins > ${proteins}.tmp
    mv ${proteins}.tmp $proteins

    # Create partitions
    /usr/local/bin/EvmUtils/partition_EVM_inputs.pl --genome $genome \\
       --gene_predictions $genes \\
       --segmentSize 2000000 \\
       --overlapSize 200000 \\
       --partition_listing $partitions \\
       --partition_dir \$PWD \\
       $protein_options \\
       $transcript_options

    # Create commands per partitions

    /usr/local/bin/EvmUtils/write_EVM_commands.pl --genome $genome \\
       --weights ${weights} \\
       --gene_predictions $genes \\
       --output_file_name $evm_out \\
       --partitions $partitions $transcript_options $protein_options > $evm_commands

    # Run commands via xargs
    cat $evm_commands \\
    | xargs \\
        -P${task.cpus} \\
        -n${task.cpus} \\
        -d'\\n' \\
        -I{} \\
        /bin/bash -c '{}'

    # Collect outputs
    /usr/local/bin/EvmUtils/recombine_EVM_partial_outputs.pl \\
        --partitions $partitions \\
        --output_file_name ${evm_out}

    # Create gffs
    /usr/local/bin/EvmUtils/convert_EVM_outputs_to_GFF3.pl \\
        --partitions $partitions \\
        --genome $genome \\
        --output_file_name ${evm_out}

    # Collect gff files
    find . \\
    | grep evm.out.gff3 \\
    | sort \\
    | xargs cat > ${meta}_evm.full.gff

    # Collect evm.out files
    find ./* \\
    | grep evm.out\$ \\
    | sort \\
    | xargs cat > ${meta}_evm.out

    # Create protein fasta
    /usr/local/bin/EvmUtils/gff3_file_to_proteins.pl ${meta}_evm.full.gff $genome prot > ${meta}_evm.full.pep

    # Create nucleotide fasta
    /usr/local/bin/EvmUtils/gff3_file_to_proteins.pl ${meta}_evm.full.gff $genome CDS > ${meta}_evm.full.cds

    # Create bed file
    /usr/local/bin/EvmUtils/gene_gff3_to_bed.pl ${meta}_evm.full.gff \\
    | sort -k1,1 -k2,2g -k3,3g > ${meta}_evm.full.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        evidencemodeler: 2.1.0
    END_VERSIONS
    """
}
