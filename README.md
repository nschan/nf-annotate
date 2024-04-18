# nf-arannotate

The current recommended workflow for assembly and annotation of _Arabidopsis_ from ONT reads is:

  * Assembly: [`nf-arassembly`](https://gitlab.lrz.de/beckerlab/nf-arassembly)
  * Annotation: This pipeline.

This pipeline is designed to annotate outputs from [`nf-arassembly`](https://gitlab.lrz.de/beckerlab/nf-arassembly).
It takes a samplesheet of genome assemblies, intitial annoations (liftoff) and *cDNA* ONT Nanopore reads.

# Graph

General Graph

```mermaid
graph TD;
    subgraph Prepare Genome
    gfasta[Genome Fasta] --> lfilt[Length filter];
    gfasta --> pseqs[Protein sequences];
    ggff[Genome GFF] --> pseqs;
    end
    subgraph Ab initio annotation
    lfilt --> AUGUSTUS;
    lfilt --> SNAP;
    lfilt --> MINIPROT;
    end
    subgraph R-Gene prediction
    pseqs --> HRP
    end
    subgraph Transcript discovery
    cDNA[cDNA Fastq] --> Porechop;
    Porechop --> minimap2;
    gfasta --> minimap2;
    ggff --> batrans[bambu transcripts]
    minimap2 --> batrans;
    batrans --> pasa;
    end
    subgraph Annotation consolidation
    AUGUSTUS --> EvidenceModeler{EvidenceModeler};
    SNAP --> EvidenceModeler;
    MINIPROT --> EvidenceModeler;
    HRP --> EvidenceModeler;
    pasa --> EvidenceModeler;
    end
    subgraph Gene Counts
    minimap2 --> bacounts[bambu counts];
    EvidenceModeler --> bacounts;
    end
    bacounts --> tsv_count>Gene Count TSV];
    subgraph Functional annotation
    EvidenceModeler --> BLASTp;
    EvidenceModeler --> pfam[Interproscan Pfam];
    BLASTp --> func[Functional annotation];
    pfam --> func;
    end
    func --> gff_anno>Annotation GFF];
    subgraph R-Gene extraction
    pfam --> rgene[R-Gene filter];
    end
    rgene --> r_tsv>R-Gene TSV];
```

Graph for HRP

```mermaid
graph TD;
  fasta[Genome Fasta] --> protseqs[Protein Sequences]
  ingff[Genome GFF] --> protseqs[Protein Sequences]
  protseqs --> pfam[Interproscan Pfam]
  pfam --> nbarc[NB-LRR extraction]
  nbarc --> meme[MEME]
  meme --> mast[MAST]
  mast --> superfam[Interproscan Superfamily]
  pfam --> rgdomains[R-Gene Identification based on Domains]
  superfam --> rgdomains
  rgdomains --> miniprot[Discovery based on known R-Genes miniprot]
  miniprot --> seqs[R-Gene sequences]
  miniprot --> rgff[R-Gene gff]
  ingff --> mergegff[Merged GFF]
  rgff --> mergegff[Merged GFF]
```

# Tubemap

![Tubemap](nf-arannotate.tubes.png)

# Usage

To run the pipeline with a samplesheet on biohpc_gen with charliecloud:

```
git clone https://gitlab.lrz.de/beckerlab/nf-arannotate
nextflow run nf-evmodeler --samplesheet 'path/to/sample_sheet.csv' \
                          --out './results' \
                          -profile biohpc_gen
```

# Parameters

| Parameter | Effect |
| --- | --- |
| `--porechop` | Run porechop on the reads? (default: `false`) |
| `--exclude_patterm` | Exclusion pattern for chromosome names (HRP, default `ATMG`, ignores mitochondrial genome) |
| `--samplesheet` | Path to samplesheet |
| `--reference_name` | Reference name (for BLAST), default: `Col-CEN` |
| `--reference_proteins` | Protein reference (defaults to Col-CEN); see known issues / blast below for additional information |
| `--r_genes` | Run R-Gene prediction pipeline?, default: `true` |
| `--min_contig_length` | minimum length of contigs to keep, default: 5000 |
| `--out` | Results directory, default: `'./results'` |

# Samplesheet

Samplesheet `.csv` with header:

```
sample,genome_assembly,liftoff,reads
```

| Column | Content |
| --- | --- |
| `sample` | Name of the sample |
| `genome_assembly` | Path to assembly fasta file |
| `liftoff` | Path to liftoff annotations |
| `reads` | Path to file containing cDNA reads |

# Procedure

This pipeline will run the following subworkflows:
  
  * `SUBSET_GENOMES`: Subset to genome to `params.min_contig_length`
  * `SUBSET_ANNOTATIONS`: Subset input gff to contigs larger than `params.min_contig_length`
  * `HRP`: Run the homology based R-gene prediction
  * `AB_INITIO`: Perform ab initio predictions:
    - `SNAP` https://github.com/KorfLab/SNAP/tree/master
    - `AUGUSTUS` https://github.com/Gaius-Augustus/Augustus (kind of paralellized)
    - `MINIPROT` https://github.com/lh3/miniprot
  * `BAMBU`: Run `porechop` (optional) on cDNA reads and align via `minimap2` in `splice:hq` mode. Then run `bambu`
  * `PASA`: Run the [PASA pipeline](https://github.com/PASApipeline/PASApipeline/wiki) on bambu output . This step starts by converting the bambu output (.gtf) by passing it through `agat_sp_convert_gxf2gxf.pl`. Subsequently transcripts are extracted (step `PASA:AGAT_EXTRACT_TRANSCRIPTS`). After running `PASApipeline` the coding regions are extracted via `transdecoder` as bundeld with pasa (`pasa_asmbls_to_training_set.dbi`)
  * `EVIDENCE_MODELER`: Take all outputs from above and the initial annotation (typically via `liftoff`) and run them through [Evidence Modeler](https://github.com/EVidenceModeler/EVidenceModeler/wiki). The implementation of this was kind of tricky, it is currently parallelized in chunks via `xargs -n${task.cpus} -P${task.cpus}`. I assume that this is still faster than running it fully sequentially. This produces the final annotations, `FUNCTIONAL` only extends this with extra information in column 9 of the gff file.
  * `GET_R_GENES`: R-Genes (NLRs) are identified in the final annotations based on `interproscan`.
  * `FUNCTIONAL`: Create functional annotations based on `BLAST` against reference and `interproscan-pfam`. Produces protein fasta. Creates `.gff` and `.gtf` outputs. Also quantifies transcripts via `bambu`.
  * `TRANSPOSONS`: Annotate transposons using `EDTA`

The weights for EVidenceModeler are defined in `assets/weights.tsv`

# Outputs

The outputs will be put into `params.out`, defaulting to `./results`. Inside the results folder, the outputs are structured according to the different subworkflows of the pipeline (`workflow/subworkflow/process`). 
All processess will emit their outputs to results.
[`AGAT`](https://github.com/NBISweden/AGAT/) is used throughout this pipeline, hopefully ensuring consistent gff formating.

# Pipeline information 

This pipeline performs a number of steps specifically aimed at discovery and annotation of NLR genes.


# Known issues & edge case handling

## Interproscan

Interproscan is run from the interproscan docker image.
The data needs to be downloaded separately and mounted into /opt/interproscan/data (see biohpc_gen.config, https://hub.docker.com/r/interpro/interproscan).
After downloading a new data-release, the container should be run once interactively to index the modles (https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html#index-hmm-models):

```bash
python3 setup.py interproscan.properties
```

## genblastG

`genblastG` was used in the original HRP publication. `genblastG` produces too many errors to be reasonably used for production tools, `miniprot` is replacing `genblastG` in this pipeline.

## BLAST / AGAT_FUNCTIONAL_ANNOTATION

`agat_sp_manage_functional_annotation.pl` is looking for `GN=` in the headers of the `.fasta` file used as a db for `BLASTP` to assign a **g**ene **n**ame.

Currently, this is handled using `sed` for a very specific case: the annotations that come with [Col-CEN-v1.2](https://github.com/schatzlab/Col-CEN).

The easiest solution would be to correctly prepare the protein fasta in such a way that it contains `GN=` with the appropriate gene names. In that case modules `MAKEBLASTDB` and `AGAT_FUNCTIONAL_ANNOTATION` need to be edited.