# . pipeline parameters



## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |

## Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `samplesheet` |  | `string` | None | True |  |
| `enable_conda` | use conda (not supported) | `boolean` |  |  |  |
| `min_contig_length` | Minimum contig length to keep | `string` | 5000 |  |  |
| `reference_proteins` | Reference proteins for miniprot | `string` | '/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/reference_genomes/Arabidopsis/Col-CEN/Col-CEN_v1.2_proteins.fasta' |  |  |
| `out` | Directory for outputs | `string` | ./results |  |  |
| `porechop` | Run porechop | `boolean` |  |  |  |
| `exclude_pattern` | Pattern to exclude for HRP | `string` | ATMG |  |  |
| `nevm` | not used | `string` | 10 |  | True |
