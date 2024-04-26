/*
==============================================
                HRP modules
==============================================
*/

include { AGAT_FILTER_BY_LENGTH } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_PROTEINS } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_MINIPROT_NLR } from '../modules/HRP/agat/main'
include { AGAT_EXTRACT_NLR } from '../modules/HRP/agat/main'
include { AGAT_COMPLEMENT } from '../modules/HRP/agat/main'
include { MEME } from '../modules/HRP/memesuite/main'
include { MAST } from '../modules/HRP/memesuite/main'
include { BEDTOOLS_GETFASTA } from '../modules/HRP/bedtools/main'
include { BEDTOOLS_CLUSTER } from '../modules/HRP/bedtools/main'
include { BEDTOOLS_NR_CLUSTERS } from '../modules/HRP/bedtools/main'
include { INTERPROSCAN as INTERPROSCAN_FULL } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_EXTENDED } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_PFAM } from '../modules/HRP/interproscan/main'
include { INTERPROSCAN_SUPERFAMILY } from '../modules/HRP/interproscan/main'
include { FILTER_R_GENES as HRP_FILTER_R_GENES } from '../modules/HRP/local/main'
include { FILTER_R_GENES_SINGLE_INPUT as FIND_R_GENES } from '../modules/HRP/local/main'
include { SEQTK_SUBSET_RPS } from '../modules/HRP/seqtk/main'
include { SEQTK_SUBSET_FL } from '../modules/HRP/seqtk/main'
include { SEQTK_SUBSET_CANDIDATES } from '../modules/HRP/seqtk/main'
include { MINIPROT_HRP } from '../modules/HRP/miniprot/main'
include { SEQKIT_GET_LENGTH } from '../modules/HRP/seqkit/main'
include { GET_R_GENE_GFF } from '../modules/HRP/local/main'

/*
==============================================
De novo annotation and EvidenceModeler modules
==============================================
*/

/*
Processing of reads, genomes and annotations
*/

include { PORECHOP } from '../modules/porechop/main.nf'
include { SEQKIT_GET_LENGTH as SEQKIT_CONTIG_LENGTH } from '../modules/seqkit/main.nf'
include { SEQTK_SUBSET_FASTA } from '../modules/seqtk/main.nf'
include { SUBSET_ANNOTATIONS } from '../modules/seqtk/main.nf'
include { AGAT_FIX_EXTRACT_TRANSCRIPTS as AGAT_EXTRACT_TRANSCRIPTS } from '../modules/agat/main.nf'
include { AGAT_GTF2GFF } from '../modules/agat/main.nf'
include { AGAT_GXF2GFF } from '../modules/agat/main.nf'
include { AGAT_GFF2GTF } from '../modules/agat/main.nf'
include { AGAT_FUNCTIONAL_ANNOTATION } from '../modules/agat/main.nf'

/*
Ab initio
*/
include { SNAP } from '../modules/snap/main.nf'
include { AUGUSTUS_PARALLEL as AUGUSTUS } from '../modules/augustus/main.nf'
include { MINIPROT } from '../modules/miniprot/main.nf'

/* 
Bambu
*/

include { ALIGN_TO_BAM as ALIGN} from '../modules/align/main.nf'
include { BAMBU } from '../modules/bambu/main'

/*
Trinity
*/
include { TRIMGALORE } from '../modules/trimgalore/main.nf'
include { ALIGN as MAP } from '../modules/star/main.nf'
include { GENOMEGENERATE } from '../modules/star/main.nf'
include { TRINITY } from '../modules/trinity/main.nf'
/*
PASA
*/
include { PASA_SEQCLEAN } from '../modules/pasa_pipeline/main.nf'
include { PASA_TD } from '../modules/pasa_pipeline/main.nf'
include { PASA_PIPELINE } from '../modules/pasa_pipeline/main.nf'
include { TRANSDECODER } from '../modules/transdecoder/main.nf'

/*
Evidence modeler
*/
include { EVIDENCEMODELER_PART_EXEC_MERGE as EV_RUN } from '../modules/evidencemodeler/main.nf'

/*
BLAST
*/

include { MAKEBLASTDB } from '../modules/blast/makeblastdb/main.nf'
include { BLASTP } from '../modules/blast/blastp/main.nf'

/*
EDTA
*/
include { AGAT_GFF2BED } from '../modules/agat/main'
include { EDTA } from '../modules/edta/main'

/* 
 ===========================================
 ===========================================
 * SUBWORKFLOWS
 ===========================================
 ===========================================
 */

workflow HRP {
    take: 
      hrp_in // tuple val(meta), path(fasta), path(gff)

    main:
      // Step 1 Extract proteins
      /* This procedure produces a few broken proteins on TAIR11.
         This is somewhat concerning, but.. well.
      */ 
      hrp_in
        .map { row -> [row[0], row[1]] }
        .set { genome }

      hrp_in
        .map { row -> [row[0], row[2]] }
        .set { ref_gff }

      AGAT_EXTRACT_PROTEINS(hrp_in, params.exclude_pattern)

      AGAT_EXTRACT_PROTEINS
        .out
        .set { proteins }
      // Step 2 Interproscan
      // This step works with spack module interproscan/5.63-95.0
      // I could not locate a container with this version.
      // Internal container build pipeline failed, the container exceeds storage on our gitlab container storage.
      //
      // It does not seem to give proper results with biocontainers/interproscan:5.59_91.0--hec16e2b_1
      // This seems to be version related.
      // For me interproscan:5.59_91.0 with -dp did not run, subjobs failed and the result was incomplete.
      // To use interproscan without -dp the version needs to be the current release.
      // As of July 2023 this is 5.63-95-0
      // I guess there are work-arounds for this, it should work with an updated container.
      INTERPROSCAN_PFAM(proteins)
      //INTERPROSCAN_EXTENDED(proteins)
      // Step 3.1 Bedfile
      proteins
        .join(INTERPROSCAN_PFAM.out.nb_bed)
        .set { bedtools_gf_in }
      // Step 3.2 Extract
      BEDTOOLS_GETFASTA(bedtools_gf_in)
      // Step 3.3 MEME
      MEME(BEDTOOLS_GETFASTA.out)
      // Step 4 MAST
      MAST(proteins
            .join(MEME.out))
      // Step 5

      proteins
        .join(INTERPROSCAN_PFAM
                .out
                .protein_tsv
                .join(
                  MAST
                  .out
                  .mast_geneids
                  )
          )
        .set { to_subset }
                    
      SEQTK_SUBSET_RPS(to_subset)

      INTERPROSCAN_SUPERFAMILY(SEQTK_SUBSET_RPS.out)
      // Step 6
      HRP_FILTER_R_GENES(
        INTERPROSCAN_PFAM
        .out
        .protein_tsv
        .join(INTERPROSCAN_SUPERFAMILY.out)
        )
      // Step 7
      SEQTK_SUBSET_FL(proteins
                      .join(
                        HRP_FILTER_R_GENES
                        .out
                        .full_length_tsv)
                      )

      // Replacing genblast with miniprot
      genome
        .join(SEQTK_SUBSET_FL.out)
        .set { miniprot_in }
  
      MINIPROT_HRP(miniprot_in) // emits gff

      AGAT_FILTER_BY_LENGTH(MINIPROT_HRP
                            .out
                            .miniprot_nlrs)

      genome
        .join(MINIPROT_HRP.out.miniprot_nlrs)
        .set { miniprot_nlr_to_extract }
      // Create proteins from gff
      AGAT_EXTRACT_MINIPROT_NLR(miniprot_nlr_to_extract)

      // Get lengths
      SEQKIT_GET_LENGTH(AGAT_EXTRACT_MINIPROT_NLR.out.extracted_nlrs)

      // Step 8.2
      BEDTOOLS_CLUSTER(AGAT_FILTER_BY_LENGTH
                        .out
                        .filtered_bed)

      // Step 8.3
      // Step 8.4
      BEDTOOLS_NR_CLUSTERS(BEDTOOLS_CLUSTER
                            .out
                            .join(SEQKIT_GET_LENGTH.out))

      // Step 9
      //   Extract annotations of non-redundant genes
      GET_R_GENE_GFF(AGAT_FILTER_BY_LENGTH
                      .out
                      .filtered_gff
                      .join(BEDTOOLS_NR_CLUSTERS.out))

      //   Extract protein sequences
      AGAT_EXTRACT_NLR(genome
                        .join(GET_R_GENE_GFF
                                .out
                                .r_genes_merged_gff))

      //   Merge R-Gene gff and input gff
      AGAT_COMPLEMENT(ref_gff
                        .join(GET_R_GENE_GFF
                                .out
                                .r_genes_merged_gff))

      AGAT_COMPLEMENT
        .out
        .merged_gff
        .set { merged_gff }

      AGAT_COMPLEMENT
        .out
        .merged_gtf
        .set { merged_gtf }

    emit:
        merged_gff
        merged_gtf

}

workflow GET_R_GENES {
    take: 
      inputs // tuple val(meta), path(fasta), path(gff)

    main:
      AGAT_EXTRACT_PROTEINS(inputs, params.exclude_pattern)

      INTERPROSCAN_PFAM(AGAT_EXTRACT_PROTEINS.out)
      INTERPROSCAN_PFAM
        .out
        .protein_tsv
        .set { pfam_out }
      
      if(params.r_genes) {
        FIND_R_GENES(INTERPROSCAN_PFAM.out.protein_tsv)
        FIND_R_GENES
          .out
          .full_length_tsv
          .set { r_gene_gff }
      }
      
    emit: 
      pfam_out
}

 workflow PREPARE_GENOMES {
  take: 
    samples // meta, fasta

  main: 
    SEQKIT_CONTIG_LENGTH(samples, params.min_contig_length)

    SEQKIT_CONTIG_LENGTH
      .out
      .contig_list
      .set { contig_lengths }

    SEQTK_SUBSET_FASTA(SEQKIT_CONTIG_LENGTH.out.large_contigs)

    SEQTK_SUBSET_FASTA
      .out
      .subset
      .set { prepared_genomes }

  emit: 
    prepared_genomes //meta, fasta
    contig_lengths // meta, lengths
 }

  workflow PREPARE_ANNOTATIONS {
  take: 
    annotations // meta, liftfoff, contig_lengths
  main: 
    SUBSET_ANNOTATIONS(annotations)

    SUBSET_ANNOTATIONS
      .out
      .annotations
      .set { annotation_subset }

  emit: 
    annotation_subset
 }
/*
 ===========================================
        RUN BAMBU
 ===========================================
*/

 workflow RUN_BAMBU {

  take: ch_bambu // sample, genome_assembly, liftoff, reads 

  main:
  
    ch_bambu
      .map { it -> [it[0], it[3]] }
      .set { ch_reads }

    ch_bambu
      .map { it -> [it[0], it[1]] }
      .set { ch_genome }

    ch_bambu
      .map { it -> [it[0], it[1], it[2]] }
      .set { ch_bambu_in }

    if(params.porechop) {

      PORECHOP(ch_reads)

      PORECHOP
        .out
        .reads
        .map { it -> [it[0],it[1]] }
        .join(ch_genome)
        .set { ch_aln }

    } else {

      ch_reads
        .join(ch_genome)
        .set { ch_aln }

    }

    ALIGN(ch_aln)
    ALIGN
      .out
      .set { alignment }
      
    BAMBU(ch_bambu_in
            .join(alignment))

    BAMBU
      .out
      .extended_gtf
      .set { bambu_gtf }

  emit:
    bambu_gtf
    alignment

 }

 /* 
 ===========================================
        AB INITIO PREDICTION
 ===========================================
 */



 workflow AB_INITIO {
  take: ch_genomes // meta, genome

  main:
    AUGUSTUS(ch_genomes)

    AUGUSTUS
      .out
      .set { augustus }

    SNAP(ch_genomes)

    SNAP
      .out
      .snap_gff
      .set { snap }

    MINIPROT(ch_genomes)
    
    MINIPROT
      .out
      .set { miniprot }

  emit:
    augustus
    snap
    miniprot
 }

 /*
 ===========================================
        TRINITY
 ===========================================
 */ 
  /*
  Add support for TRINITY, to enable annotation with short-reads
  */

 /*
 Accessory function to create input for trinity
 modified from nf-core/rnaseq/subworkflows/local/input_check.nf
 */

  def create_shortread_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id       = row[0]
    meta.paired   = row.paired.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def shortreads = []
    if (!file(row.shortread_F).exists()) {
        exit 1, "ERROR: shortread_F fastq file does not exist!\n${row.shortread_F}"
    }
    if (!meta.paired) {
        shortreads = [ meta.id, meta.paired, [ file(row.shortread_F) ] ]
    } else {
        if (!file(row.shortread_R).exists()) {
            exit 1, "ERROR: shortread_R fastq file does not exist!\n${row.shortread_R}"
        }
        shortreads = [ meta.id, meta.paired, [ file(row.shortread_F), file(row.shortread_R) ] ]
    }
    return shortreads
  }

 workflow STAR {
  take:
    star_in // sample, genome, gff, paired, reads]
  
  main:
    star_in
     .map { it -> [it[0], it[2]] }
     .set { gff_file }
    star_in
      .map { it -> [it[0], it[3], it[4]] }
      .set { reads } // sample, paired, reads
    AGAT_GFF2GTF(gff_file)
    star_in
      .map { it -> [it[0], it[1]] }
      .join(AGAT_GFF2GTF.out)
      .set { star_genome }

    GENOMEGENERATE(star_genome) // sample, index

    GENOMEGENERATE
      .out
      .join(AGAT_GFF2GTF.out) // sample, index, gtf
      .join(reads) // sample, index, gtf, paired, reads
      .set { map_in }
    MAP(map_in)
    MAP.out.bam_sorted.set { sorted_bam }
  
  emit:
    sorted_bam
 }

 workflow RUN_TRINITY {
  take:
    ch_input // define samplesheet: sample,genome, gff, short1, short2, paired
  
  main:
    ch_input
      .map(create_shortread_channel)
      .set { ch_short_reads }
    
    // Trim, takes: [sample, paired, [reads]]
    TRIMGALORE(ch_short_reads)
    TRIMGALORE
      .out
      .reads
      .set { trimmed_reads }
    ch_input
      .map { it -> [ it[0], it[1], it[2] ] } // sample, genome, gff
      .join( trimmed_reads ) // sample, genome, gff, paired, reads
      .set { splice_aln_in } 
    // STAR, takes: [sample, genome, gff, paired, reads]
    STAR(splice_aln_in)
    // Tinity, takes: [sample, bam]
    TRINITY(STAR.out)
    TRINITY.out.transcript_fasta.set{ transcripts }

  emit:
    transcripts
 }
 /*
 ===========================================
        PASA
 ===========================================
 */
 workflow PASA {
  take: 
    ch_genomes
    ch_transcripts

  main:
  // trinity gives fasta file
  if(params.short_reads) {
    ch_transcripts
      .set { transcripts }
  }

  // bambu gives a gtf file
  if(!params.short_reads) {
    AGAT_GTF2GFF(ch_transcript_gff)
    AGAT_EXTRACT_TRANSCRIPTS(AGAT_GTF2GFF.out)
    AGAT_EXTRACT_TRANSCRIPTS.out.set{ transcripts }
  }

    PASA_PIPELINE(ch_genomes.join(transcripts))

    PASA_PIPELINE
      .out
      .pasa_assembly_fasta
      .join(PASA_PIPELINE
              .out
              .pasa_assembly_gff)
      .set { td_in }
    
    PASA_PIPELINE
      .out
      .pasa_assembly_gff
      .set { pasa }

    PASA_TD(td_in)

    PASA_TD
      .out
      .genome_gff
      .set { pasa_td }

  emit:
    pasa
    pasa_td
 }
 /* 
 ===========================================
        CDS from liftoff
 ===========================================
 */
workflow CDS_FROM_ANNOT {
    take: ch_annot

    main: 
    AGAT_EXTRACT_TRANSCRIPTS(ch_annot)

    TRANSDECODER(AGAT_EXTRACT_TRANSCRIPTS.out)

    TRANSDECODER
      .out
      .td_gff
      .set { transdecoder }

    emit:
    transdecoder
}

/*
This workflow was adapted from nf-core/genomeannotator
The simple EVIDENCEMODELER module I created fails for unclear
reasons.
*/

workflow EV_MODELER {
    take:
      ev_in

    main:
    EV_RUN(ev_in)

    emit: 
    EV_RUN.out.gff
}

 /* 
 ===========================================
        BLAST against references
 ===========================================
 */

 workflow BLAST {
    take:
      query // meta, genome, annotations
      reference // meta, protein fasta

    main:
      MAKEBLASTDB(reference)
      AGAT_EXTRACT_PROTEINS(query, params.exclude_pattern)
      BLASTP(AGAT_EXTRACT_PROTEINS.out,
             MAKEBLASTDB.out.db.first(),
             "tsv" )
      BLASTP
        .out
        .tsv
        .set { blast_table }

    emit:
      blast_table
 } 

  /* 
 ===========================================
      Combine functional annotations
 ===========================================
 */

 workflow FUNCTIONAL {
  take:
    annotations // meta, gff
    blast_output // meta, tsv
    interpro_tsv // meta, tsv
    blast_reference // meta2, protein fasta
    genomes // meta, fasta
    alignments // meta, bam
  
  main:
    annotations
      .join(blast_output)
      .join(interpro_tsv)
      .set { annotation_and_function }
    
    AGAT_FUNCTIONAL_ANNOTATION(annotation_and_function, blast_reference)

    AGAT_GFF2GTF(AGAT_FUNCTIONAL_ANNOTATION.out.gff_file)
 
    if(!params.short_reads) {
      genomes
        .join(AGAT_GFF2GTF.out)
        .join(alignments)
        .set { bambu_in }
      BAMBU(bambu_in) //takes: tuple val(meta), path(fasta), path(gtf), path(bams)
    }
    

 }

  /* 
 ===========================================
        Annotate transpososons
 ===========================================
 */

 workflow TRANPOSONS {
  take:
    annotated_genome // meta, fasta, gff
      
  main:
    
    annotated_genome
      .map { it -> [ it[0], it[1]] }
      .set { genome }

    annotated_genome
      .map { it -> [ it[0], it[2]] }
      .set { annotations }
    
    AGAT_GFF2BED(annotations)

    AGAT_EXTRACT_TRANSCRIPTS(annotated_genome)

    genome
      .join(AGAT_GFF2BED.out)
      .join(AGAT_EXTRACT_TRANSCRIPTS.out)
      .set { edta_in }

    EDTA(edta_in)

    emit:
      EDTA.out.transposon_annotations
      EDTA.out.transposon_summary
 }