# nf-arannotate

Annotate outputs from [`nf-arassembly`](https://gitlab.lrz.de/beckerlab/nf-arassembly).
It takes a samplesheet of genome assemblies, intitial annoations (liftoff) and cDNA reads.
This pipeline is a combination of [`nf-hrp`](https://gitlab.lrz.de/beckerlab/nf-hrp) and [`nf-evmodeler`](https://gitlab.lrz.de/beckerlab/nf-evmodeler)


# Usage

To run the pipeline with a samplesheet on biohpc_gen:

```
git clone https://gitlab.lrz.de/beckerlab/nf-arannotate
nextflow run nf-evmodeler --samplesheet 'path/to/sample_sheet.csv' \
                          --out './results' \
                          -profile charliecloud,biohpc_gen
```

# Parameters

| Parameter | Effect |
| --- | --- |
| `--porechop` | Run porechop on the reads? (default: `false`) |
| `--exclude_patterm` | Exclusion pattern for chromosome names (HRP, default `ATMG`, ignores mitochondrial genome) |
| `--samplesheet` | Path to samplesheet |
| `--protein_ref` | Protein reference (defaults to Col-CEN) |
| `--min_contig_length` | minimum length of contigs to keep (default: 5000) |
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

This pipeline will:
  
  * `SUBSET_GENOMES`: Subset to genome to `params.min_contig_length`
  * `SUBSET_ANNOTATIONS`: Subset input gff to contigs larger than `params.min_contig_length`
  * `HRP`: Run the homology based R-gene prediction
  * `AB_INITIO`: Perform ab initio predictions:
    - `SNAP` https://github.com/KorfLab/SNAP/tree/master
    - `AUGUSTUS` https://github.com/Gaius-Augustus/Augustus (kind of paralellized)
    - `MINIPROT` https://github.com/lh3/miniprot
  * `BAMBU`: Run porechop on cDNA reads and align via minimap2. Then run `bambu`
  * `PASA`: Run the PASA pipeline on bambu output (https://github.com/PASApipeline/PASApipeline/wiki). This step starts by converting the bambu output (.gtf) by passing it through `agat_sp_convert_gxf2gxf.pl`. Subsequently transcripts are extracted (step `PASA:AGAT_EXTRACT_TRANSCRIPTS`). After running `PASApipeline` the coding regions are extracted via `transdecoder` as bundeld with pasa (`pasa_asmbls_to_training_set.dbi`)
  * `EVIDENCE_MODELER`: Take all outputs from above and the initial annotation (typically via `liftoff`) and run them through Evidence Modeler (https://github.com/EVidenceModeler/EVidenceModeler/wiki). The implementation of this was kind of difficult, it is currently parallelized in chunks via `xargs -n${task.cpus} -P${task.cpus}`. I assume that this is still faster than running it fully sequentially.

The weights for EVidenceModeler are defined in `assets/weights.tsv`

# Outputs

The outputs will be put into `params.out`, defaulting to `./results`. Inside the results folder, the outputs are structured according to the different subworkflows of the pipeline (`workflow/subworkflow/process`). 
All processess will emit their outputs to results, for the ab initio predictions those are gff files, also for pasa. Evidence modeler will emit gff, bed, pep and cds files containing annotations (gff & bed) and translated sequences (.pep and .cds are fasta formatted).

# Additional information

The current recommended workflow for assembly and annotation of Arabidopsis from ONT reads is:

  * Assembly: `nf-arassembly` [https://gitlab.lrz.de/beckerlab/nf-arassembly]
  * Annotation: This pipeliine.