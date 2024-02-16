

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRANSDECODER {
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(meta), path(transcripts)

    output:
        tuple val(meta), path("${transcripts}.transdecoder_dir/*.pep"), emit: td_pep
        tuple val(meta), path("${transcripts}.transdecoder_dir/*.gff3"), emit: td_gff
        tuple val(meta), path("${transcripts}.transdecoder_dir/*.cds"), emit: td_cds      

    script:
        def workflow="${task.process}".tokenize(":")[-2].toLowerCase()
        """
        TransDecoder.LongOrfs -t ${transcripts} 
        mv ${transcripts}.transdecoder_dir/longest_orfs.gff3 ${transcripts}.transdecoder_dir/${meta}_${workflow}_longest_orfs.gff3 
        mv ${transcripts}.transdecoder_dir/longest_orfs.pep ${transcripts}.transdecoder_dir/${meta}_${workflow}_longest_orfs.pep
        mv ${transcripts}.transdecoder_dir/longest_orfs.cds ${transcripts}.transdecoder_dir/${meta}_${workflow}_longest_orfs.cds 
        """
}