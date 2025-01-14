[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12759772.svg)](https://zenodo.org/doi/10.5281/zenodo.12759772)

The goal of [`genomeassembler`](https://github.com/nschan/genomeassembler) and [`nf-annotate`](https://github.com/nschan/nf-annotate) is to make to genome assembly and annotation workflows accessible for a broader community, particularily for plant-sciences. Long-read sequencing technologies are already cheap and will continue to drop in price, genome sequencing will soon be available to many researchers without a strong bioinformatic background. 
The assembly is naturally quite organisms agnostic, but the annotation pipeline contains some steps that may not make sense for other eukaryotes, unless there is a particular interest in NB-LRR genes.

# nf-arannotate

The current recommended workflow for assembly and annotation of _Arabidopsis_ from long reads is:

  * Assembly: [`genomeassembler`](https://github.com/nschan/genomeassembler)
  * Annotation: This pipeline.

This pipeline is designed to annotate outputs from [`nf-genomeassembly`](https://github.com/nschan/nf-genomeassembly).
It takes a samplesheet of genome assemblies, intitial annotations (liftoff) and *cDNA* ONT Nanopore reads or pacbio isoseq reads. If no long transcriptome reads are available short reads can also be used.

If `--short_reads` is `true` the pipeline takes short reads instead of long cDNA. This is probably better than no reads, but for high-quality annotations long transcriptome reads are recommended.

# Usage

To run the pipeline with a samplesheet on biohpc_gen with charliecloud:

```
git clone https://github.com/nschan/nf-annotate
nextflow run nf-annotate --samplesheet 'path/to/sample_sheet.csv' \
                          --out './results' \
                          -profile biohpc_gen
```

# Parameters

| Parameter | Effect |
| --- | --- |
| `--samplesheet` | Path to samplesheet |
| `--preprocess_reads` | Run `porechop` on ONT reads or [`LIMA`-`REFINE`](https://isoseq.how/getting-started.html) on pacbio reads? (default: `false`) |
| `--exclude_pattern` | Exclusion pattern for chromosome names (HRP, default `ATMG`, ignores mitochondrial genome) |
| `--reference_name` | Reference name (for BLAST), default: `Col-CEN` |
| `--reference_proteins` | Protein reference (defaults to Col-CEN); see known issues / blast below for additional information |
| `--gene_id_pattern` | Regex to capture gene name in initial annoations. Default: ` "AT[1-5C]G[0-9]+.[0-9]+|evm[0-9a-z\\.]*|ATAN.*" ` will capture TAIR IDs, evm IDs and ATAN  |
| `--r_genes` | Run R-Gene prediction pipeline?, default: `true` |
| `--augustus_species` | Species to for agustus, default: `"arabidopsis"` |
| `--snap_organism` | Model to use for snap, default: `"A.thaliana"` |
| `--mode` | Specify `'ont'` or `'pacbio'`. Default `'ont'` |
| `--aligner` | Aligner for long-reads. Options are `'minimap2'` or `ultra`. Default: `'minimap2'` |
| `--pacbio_polya` | Require (and trim) polyA tails from pacbio reads? Default: `true` |
| `--primers` | File containing primers used for pacbio sequencing (required if `--mode` is 'pacbio'). Default : `null` |
| `--short_reads` | Provide this parametere if the transcriptome reads are short reads (see below). Default: `false` |
| `--bamsortram` | *Short-reads only*: passed to STAR for `--limitBAMsortRAM`. Specifies RAM available for BAM sorting, in bytes. Default: `0` |
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

If `--short_reads` is used the samplesheet should look like:

```
sample,genome_assembly,liftoff,paired,shortread_F,shortread_R
sampleName,assembly.fasta,reference.gff,true,short_F1.fastq,short_F2.fastq
```

| Column | Content |
| --- | --- |
| `sample` | Name of the sample |
| `genome_assembly` | Path to assembly fasta file |
| `liftoff` | Path to liftoff annotations |
| `pair` | `true` or `false` depending on whether the short reads are paired |
| `shortread_F` | Path to forward reads |
| `shortread_R` | Path to reverse reads |

> If there is only one type of read shortread_R should remain empty and paired should be `false`

> NB: It is possible to mix paired and unpaired reads within one samplesheet, e.g. when performing annotation of many genomes with heterogenious data availability.

> NB: It is *not* possible to mix long and short reads in a single samplesheet.

# Procedure

This pipeline will run the following subworkflows:
  
  * `SUBSET_GENOMES`: Subset to genome to `params.min_contig_length`
  * `SUBSET_ANNOTATIONS`: Subset input gff to contigs larger than `params.min_contig_length`
  * `HRP`: Run the homology based R-gene prediction
  * `AB_INITIO`: Perform ab initio predictions:
    - `SNAP` https://github.com/KorfLab/SNAP/tree/master
    - `AUGUSTUS` https://github.com/Gaius-Augustus/Augustus (kind of paralellized)
    - `MINIPROT` https://github.com/lh3/miniprot
  * `BAMBU` (long cDNA reads): Run `porechop` (optional) on cDNA reads. These reads are aligned via `minimap2` in `splice:hq` mode or using `ultra`, depending on the value of `params.aligner`. Then run `bambu`
  * `TRINITY` (short cDNA reads): Run `Trim Galore!` on the short reads, followed by `STAR` for alignment and `TRINITY` for transcript discovery from the alignment.
  * `PASA`: Run the [PASA pipeline](https://github.com/PASApipeline/PASApipeline/wiki) on bambu output . This step starts by converting the bambu output (.gtf) by passing it through `agat_sp_convert_gxf2gxf.pl`. Subsequently transcripts are extracted (step `PASA:AGAT_EXTRACT_TRANSCRIPTS`). After running `PASApipeline` the coding regions are extracted via `transdecoder` as bundeld with pasa (`pasa_asmbls_to_training_set.dbi`)
  * `EVIDENCE_MODELER`: Take all outputs from above and the initial annotation (typically via `liftoff`) and run them through [Evidence Modeler](https://github.com/EVidenceModeler/EVidenceModeler/wiki). The implementation of this was kind of tricky, it is currently parallelized in chunks via `xargs -n${task.cpus} -P${task.cpus}`. I assume that this is still faster than running it fully sequentially. This produces the final annotations, `FUNCTIONAL` only extends this with extra information in column 9 of the gff file.
  * `GET_R_GENES`: R-Genes (NLRs) are identified in the final annotations based on `interproscan`.
  * `FUNCTIONAL`: Create functional annotations based on `BLAST` against reference and `interproscan-pfam`. Produces protein fasta. Creates `.gff` and `.gtf` outputs. Also quantifies transcripts via `bambu`.
  * `TRANSPOSONS`: Annotate transposons using `HiTE`

The weights for EVidenceModeler are defined in `assets/weights.tsv`

# Outputs

The outputs will be put into `params.out`, defaulting to `./results`. Inside the results folder, the outputs are structured according to the different subworkflows of the pipeline (`workflow/subworkflow/process`). 
All processess will emit their outputs to results.
[`AGAT`](https://github.com/NBISweden/AGAT/) is used throughout this pipeline, hopefully ensuring consistent gff formating.

# Graph

Graph for HRP

```mermaid
graph TD;
  fasta>Genome Fasta] --> protseqs[Protein Sequences]
  ingff>Genome GFF] --> protseqs[Protein Sequences]
  protseqs --> pfam[Interproscan Pfam]
  pfam --> nbarc[NB-LRR extraction]
  nbarc --> meme[MEME]
  meme --> mast[MAST]
  mast --> superfam[Interproscan Superfamily]
  pfam --> rgdomains[R-Gene Identification based on Domains]
  superfam --> rgdomains
  rgdomains --> miniprot[miniprot: discovery based on known R-genes]
  miniprot --> seqs>R-Gene sequences]
  miniprot --> rgff[R-Gene gff]
  ingff --> mergegff>Merged GFF]
  rgff --> mergegff
```

## Overall graph

```mermaid
%%{init: {'theme': 'dark',
          'themeVariables':{
            'commitLabelColor': '#cccccc',
            'commitLabelBackground': '#434244',
            'commitLabelFontSize': '12px',
            'tagLabelFontSize': '12px',
            'git0': '#8db7d2',
            'git1': '#58508d',
            'git2': '#bc5090',
            'git3': '#ff6361',
            'git4': '#ffa600',
            'git5': '#74a892',
            'git6': '#d69e49',
            'git7': '#00ffff'
            },
          'gitGraph': {
            'mainBranchName': "Prepare Genome",
             'parallelCommits': false
             } 
          }
}%%

gitGraph TB:
  commit id: "Genome fasta"
  commit id: "Length filter [seqtk]" tag: "fasta"
  branch "HRP"
  branch "Ab initio<br>prediction"
  branch "Transcript<br>discovery"
  branch "Evidence Modeler"
  checkout "Prepare Genome"
  commit id: "Protein sequences [agat]"
  checkout "HRP"
  commit id: "NLR Extraction"
  commit id: "InterproScan PFAM"
  commit id: "MEME"
  commit id: "MAST"
  commit id: "InterproScan Superfamily"
  commit id: "Genome scan [miniprot]"
  commit id: "Merge with input"
  checkout "Evidence Modeler"
  merge "HRP" tag: "R-gene GFF"
  checkout "Ab initio<br>prediction"
  commit id: "AUGUSTUS"
  checkout "Evidence Modeler"
  merge "Ab initio<br>prediction" tag: "AUGUSTUS GFF"
  checkout "Ab initio<br>prediction"
  commit id: "SNAP"
  checkout "Evidence Modeler"
  merge "Ab initio<br>prediction" tag: "SNAP GFF"
  checkout "Ab initio<br>prediction"
  commit id: "miniprot"
  checkout "Evidence Modeler"
  merge "Ab initio<br>prediction" tag: "miniprot GFF"
  checkout "Transcript<br>discovery"
  commit id: "Reads" tag: "fasta"
  commit id: "Porechop / Trim Galore"
  commit id: "minimap2 / STAR"
  commit id: "bambu / Trinity"
  checkout "Evidence Modeler"
  merge "Transcript<br>discovery" tag: "Transcript GFF"
  commit type: HIGHLIGHT id: "Merged GFF"
  branch "Functional<br>annotation"
  branch "Tranposon<br>annotation"
  checkout "Functional<br>annotation"
  commit id: "BLAST"
  commit id: "InterproScan"
  commit id: "Functional annotation [agat]" tag: "Gene GFF" type: HIGHLIGHT
  checkout "Tranposon<br>annotation"
  commit type: HIGHLIGHT id: "HiTE" tag: "Transposon GFF"
```

# Tubemap

![Tubemap](nf-annotate.map.png)

# Pipeline information 

This pipeline performs a number of steps specifically aimed at discovery and annotation of NLR genes.

# Output

The pipeline will put all outputs into the directory specified via `--out` (default: `./results`)

# Known issues & edge case handling

## Interproscan

Interproscan is run from the interproscan docker image.
The data needs to be downloaded separately and mounted into /opt/interproscan/data (see biohpc_gen.config, https://hub.docker.com/r/interpro/interproscan).
After downloading a new data-release, the container should be run once interactively to index the modles (https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html#index-hmm-models):

```bash
python3 setup.py interproscan.properties
```

## BLAST / AGAT_FUNCTIONAL_ANNOTATION

`agat_sp_manage_functional_annotation.pl` is looking for `GN=` in the headers of the `.fasta` file used as a db for `BLASTP` to assign a **g**ene **n**ame.

Currently, this is handled using `sed` for a very specific case: the annotations that come with [Col-CEN-v1.2](https://github.com/schatzlab/Col-CEN).

The easiest solution would be to correctly prepare the protein fasta in such a way that it contains `GN=` with the appropriate gene names. In that case modules `MAKEBLASTDB` and `AGAT_FUNCTIONAL_ANNOTATION` need to be edited.