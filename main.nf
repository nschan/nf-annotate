nextflow.enable.dsl = 2 
params.samplesheet = false
params.enable_conda = false
params.porechop = false
params.publish_dir_mode = 'copy'
params.min_contig_length = 5000
params.exclude_pattern = "ATMG"
params.reference_name = "Col-CEN"
params.reference_proteins = '/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/reference_genomes/Arabidopsis/Col-CEN/Col-CEN_v1.2_proteins.fasta'
params.augustus_species = "arabidopsis"
params.gene_id_pattern = "AT[1-5C]G[0-9]+.[0-9]+|evm[0-9a-z\\.]*|ATAN.*"
params.r_genes = true
params.short_reads = false
params.out = './results'
params.nevm = 10

include { HRP } from './subworkflows/main.nf'
include { GET_R_GENES } from './subworkflows/main.nf'
include { PREPARE_GENOMES } from './subworkflows/main.nf'
include { PREPARE_ANNOTATIONS } from './subworkflows/main.nf'
include { AB_INITIO } from './subworkflows/main.nf'
include { RUN_BAMBU } from './subworkflows/main.nf'
include { RUN_TRINITY } from './subworkflows/main.nf'
include { PASA } from './subworkflows/main.nf'
include { CDS_FROM_ANNOT } from './subworkflows/main.nf'
include { EV_MODELER } from './subworkflows/main.nf'
include { EDTA } from './subworkflows/main.nf'
include { EDTA_ANNOTATE } from './subworkflows/main.nf'
include { TRANPOSONS } from './subworkflows/main.nf'
include { BLAST } from './subworkflows/main.nf'
include { FUNCTIONAL } from './subworkflows/main.nf'
include { AGAT_GXF2GFF } from './modules/agat/main.nf'
include { create_shortread_channel } from './subworkflows/main.nf'

log.info """\
==============================================================================================================================================
==============================================================================================================================================

███▄▄▄▄      ▄████████    ▄████████ ███▄▄▄▄   ███▄▄▄▄    ▄██████▄      ███        ▄████████     ███        ▄████████ 
███▀▀▀██▄   ███    ███   ███    ███ ███▀▀▀██▄ ███▀▀▀██▄ ███    ███ ▀█████████▄   ███    ███ ▀█████████▄   ███    ███ 
███   ███   ███    █▀    ███    ███ ███   ███ ███   ███ ███    ███    ▀███▀▀██   ███    ███    ▀███▀▀██   ███    █▀  
███   ███  ▄███▄▄▄       ███    ███ ███   ███ ███   ███ ███    ███     ███   ▀   ███    ███     ███   ▀  ▄███▄▄▄     
███   ███ ▀▀███▀▀▀     ▀███████████ ███   ███ ███   ███ ███    ███     ███     ▀███████████     ███     ▀▀███▀▀▀     
███   ███   ███          ███    ███ ███   ███ ███   ███ ███    ███     ███       ███    ███     ███       ███    █▄  
███   ███   ███          ███    ███ ███   ███ ███   ███ ███    ███     ███       ███    ███     ███       ███    ███ 
 ▀█   █▀    ███          ███    █▀   ▀█   █▀   ▀█   █▀   ▀██████▀     ▄████▀     ███    █▀     ▄████▀     ██████████ 
                                                                                                                     
                                                                                                                                                          
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
     find R genes    : ${params.r_genes}
     short reads     : ${params.short_reads}
   outdir            : ${params.out}
   conda             : ${params.enable_conda}

==============================================================================================================================================
==============================================================================================================================================
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

    /*
    New samplesheet
    sample, genome_assembly, liftoff, read_type, forward, reverse
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
    
    AB_INITIO(ch_genomes)

    if(!params.r_genes) {
      AGAT_GFF2GTF(PREPARE_ANNOTATIONS
        .out
        .annotation_subset)
      AGAT_GFF2GTF
        .out
        .set { ch_annotation_subset }
    }

    if(params.r_genes) {
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
        .merged_gtf
        .set { ch_annotation_subset }
    }

    ch_genomes
      .join(ch_annotation_subset)
      .join(ch_samples
            .map { it -> [it.sample, it.reads] } 
            ) 
      .set { ch_bambu } // sample, genome_assembly, liftoff, reads 
    
    /*
    Attempt to add branching:

    New samplesheet
    sample, genome_assembly, liftoff, paired, forward, reverse
    
    ch_samples
      .branch {
        sample, genome_assembly, liftoff, paired, forward, reverse ->
          longread: paired == "long" // everything that is not long is short
          shortread: true // short-reads can be true or false for paired
     }
      .set { ch_samples }

    // Long reads go into bambu

    ch_genomes
      .join(ch_annotation_subset)
      .join(ch_samples.longread
            .map { it -> [ it.sample, it.forward ] } 
            ) 
      .set { ch_bambu }
          RUN_BAMBU(ch_bambu)

    RUN_BAMBU
      .out
      .alignment
      .set { cdna_alignment }
    
    ch_genomes
      .join(
        RUN_BAMBU
        .out
        .bambu_gtf
        )
      .set { long_transcripts }

    // Short reads go into trinity
      RUN_TRINITY(ch_samples.shortread)
      RUN_TRINITY
        .out
        .set { short_transcripts }

    long_transcripts.mix(short_transcripts)
      .set { transcripts }
    */

    // Transcript discovery
    // long reads
    if(!params.short_reads) {
    RUN_BAMBU(ch_bambu)

    RUN_BAMBU
      .out
      .alignment
      .set { cdna_alignment }
    
    ch_genomes
      .join(
        RUN_BAMBU
        .out
        .bambu_gtf
        )
      .set { transcripts }
    }
    // short reads
    if(params.short_reads) {

      Channel.empty().set { cdna_alignment }
      // create input for trinity based on subset genomes
      ch_genomes
        .join(PREPARE_ANNOTATIONS
               .out
               .annotation_subset)
        .join(ch_samples
               .map( it -> [ it.sample, it.shortread_F, it.shortread_R, it.paired ] ))
        .set { ch_trinity }
           
      RUN_TRINITY(ch_trinity)

      RUN_TRINITY
        .out
        .set { transcripts }
    }

    
    PASA(ch_genomes, transcripts)
    
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
        .set { ch_evm_annotations } // tuple val(meta), path(fasta), path(gff)
        
    GET_R_GENES(ch_evm_annotations) 

    GET_R_GENES
      .out
      .pfam_out
      .set { ch_interproscan_results }

    Channel
      .of( [params.reference_name, params.reference_proteins] )
      .set { ch_ref_proteins }

    BLAST(ch_evm_annotations, ch_ref_proteins)

    FUNCTIONAL(
      EV_MODELER.out,
      BLAST.out.blast_table,
      ch_interproscan_results,
      ch_ref_proteins.first(),
      ch_genomes,
      cdna_alignment
    )

    //EDTA(ch_genomes)
    //EDTA_ANNOTATE(ch_evm_annotations, EDTA.out)
    TRANPOSONS(ch_evm_annotations)

 }

 workflow {
   ANNOTATE()
 }