nextflow.enable.dsl = 2 
params.samplesheet = false
params.enable_conda = false
params.porechop = false
params.publish_dir_mode = 'copy'
params.min_contig_length = 5000
params.exclude_pattern = "ATMG"
params.reference_proteins = '/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/reference_genomes/Arabidopsis/Col-CEN/Col-CEN_v1.2_proteins.fasta'
params.out = './results'
params.nevm = 10

include { HRP } from './subworkflows/main.nf'
include { GET_R_GENES } from './subworkflows/main.nf'
include { PREPARE_GENOMES } from './subworkflows/main.nf'
include { PREPARE_ANNOTATIONS } from './subworkflows/main.nf'
include { RUN_BAMBU } from './subworkflows/main.nf'
include { AB_INITIO } from './subworkflows/main.nf'
include { PASA } from './subworkflows/main.nf'
include { CDS_FROM_ANNOT } from './subworkflows/main.nf'
include { EV_MODELER } from './subworkflows/main.nf'

log.info """\
==============================================================================================================================================
==============================================================================================================================================

███▄▄▄▄      ▄████████    ▄████████    ▄████████    ▄████████ ███▄▄▄▄   ███▄▄▄▄    ▄██████▄      ███        ▄████████     ███        ▄████████ 
███▀▀▀██▄   ███    ███   ███    ███   ███    ███   ███    ███ ███▀▀▀██▄ ███▀▀▀██▄ ███    ███ ▀█████████▄   ███    ███ ▀█████████▄   ███    ███ 
███   ███   ███    █▀    ███    ███   ███    ███   ███    ███ ███   ███ ███   ███ ███    ███    ▀███▀▀██   ███    ███    ▀███▀▀██   ███    █▀  
███   ███  ▄███▄▄▄       ███    ███  ▄███▄▄▄▄██▀   ███    ███ ███   ███ ███   ███ ███    ███     ███   ▀   ███    ███     ███   ▀  ▄███▄▄▄     
███   ███ ▀▀███▀▀▀     ▀███████████ ▀▀███▀▀▀▀▀   ▀███████████ ███   ███ ███   ███ ███    ███     ███     ▀███████████     ███     ▀▀███▀▀▀     
███   ███   ███          ███    ███ ▀███████████   ███    ███ ███   ███ ███   ███ ███    ███     ███       ███    ███     ███       ███    █▄  
███   ███   ███          ███    ███   ███    ███   ███    ███ ███   ███ ███   ███ ███    ███     ███       ███    ███     ███       ███    ███ 
 ▀█   █▀    ███          ███    █▀    ███    ███   ███    █▀   ▀█   █▀   ▀█   █▀   ▀██████▀     ▄████▀     ███    █▀     ▄████▀     ██████████ 
                                      ███    ███                                                                                               
                                                                                                                                                          
----------------------------------------------------------------------------------------------------------------------------------------------
Niklas Schandry                                  niklas@bio.lmu.de                               https://gitlab.lrz.de/beckerlab/nf-arannotate                                          
----------------------------------------------------------------------------------------------------------------------------------------------

  Parameters:
     samplesheet     : ${params.samplesheet}
     min contig size : ${params.min_contig_length}
     run porechop    : ${params.porechop}
     protein ref     : 
      ${params.reference_proteins}
     exlude_pattern  : ${params.exclude_pattern}
   outdir            : ${params.out}
   conda             : ${params.enable_conda}

=====================================================================================================================================
=====================================================================================================================================
"""
    .stripIndent(false)


/* 
 ===========================================
 ===========================================
 * MAIN WORKFLOW
 ===========================================
 ===========================================
 */
 workflow ANNOTATE {
  
    if(params.enable_conda) {
      exit 1, "conda is not supported, please use containers"
    }
    if(params.samplesheet) {
    Channel.fromPath(params.samplesheet) 
           .splitCsv(header:true)
           .set { ch_samples }
    }
    else {
    exit 1, 'Input samplesheet not specified!'
    }
    /*
    Samplesheet:
    sample,genome_assembly,liftoff,reads
    */
    ch_samples
      .map { it -> [it.sample, it.genome_assembly] }
      .set { ch_genomes_subset }

    ch_samples
      .map { it -> [it.sample, it.liftoff] }
      .set { ch_initial_annotations }

    PREPARE_GENOMES(ch_genomes_subset)
    PREPARE_GENOMES
      .out
      .prepared_genomes
      .set { ch_genomes }

    PREPARE_ANNOTATIONS(ch_initial_annotations
                        .join(
                          PREPARE_GENOMES
                          .out
                          .contig_lengths
                          )
                        )
    ch_genomes
      .join(
        PREPARE_ANNOTATIONS
        .out
        .annotation_subset
           )
      .set { ch_hrp_in }

    HRP(ch_hrp_in)
   
    HRP
      .out
      .merged_gff
      .set { ch_annotation_subset }

    ch_genomes
      .join(ch_annotation_subset)
      .join(ch_samples.map { it -> [it.sample, it.reads] } ) 
      .set { ch_bambu } // sample, genome_assembly, liftoff, reads 

    AB_INITIO(ch_genomes)

    RUN_BAMBU(ch_bambu)
    
    ch_genomes
      .join(
        RUN_BAMBU
        .out
        .bambu_gtf
        )
      .set { ch_genomes_bambu }
    
    PASA(ch_genomes, ch_genomes_bambu)
    
    AGAT_GXF2GFF(
      ch_annotation_subset
      .join(AB_INITIO
            .out
            .augustus
            )
      .join(AB_INITIO
            .out
            .snap
            )
      .join(PASA
            .out
            .pasa
            )
      )
                  
    ch_genomes
        .join(AGAT_GXF2GFF
              .out
              )
        .join(AB_INITIO
              .out
              .miniprot
              )
        .join(PASA
              .out
              .pasa_td
              )
        .set { ch_em_in }

    EV_MODELER(ch_em_in)

    ch_genomes
        .join(EV_MODELER
              .out
             )
        .set { ch_evm_annotations }

    GET_R_GENES(ch_evm_annotations)

 }

 workflow {
   ANNOTATE()
 }