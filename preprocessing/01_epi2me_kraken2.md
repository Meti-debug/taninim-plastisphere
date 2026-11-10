# Taxonomic classification: EPI2ME wf-metagenomics

Classification was run in the EPI2ME desktop application, not from a shell, so
this file records the run instead of reproducing it as a script.

| setting | value |
|---|---|
| workflow | `wf-metagenomics` |
| workflow version | 2.14.2 |
| classifier | Kraken2, with Bracken re-estimation |
| database | `Standard-8` (RefSeq archaea, bacteria, viruses, plasmids, UniVec_Core) |
| `taxonomic_rank` | S (species) |
| `kraken2_confidence` | 0 |
| `min_len` | 0 |
| `abundance_threshold` | 0 |
| `bracken_threshold` | 10 |
| input | `fastq_pass`, one directory per barcode |

All four threshold parameters are the workflow defaults. No additional length
or quality filtering was applied before classification: the reads were the
Q >= 9 pass reads written by the basecaller.

Read concatenation per barcode inside the workflow is done by fastcat v0.24.2.

The workflow exports a species-level count table,
`wf-metagenomics-counts-species.csv`. Convert it with:

```bash
python preprocessing/02_epi2me_counts_to_table.py \
    wf-metagenomics-counts-species.csv barcode_map.tsv
```

`barcode_map.tsv` is two tab-separated columns, barcode and library name.

To repeat the classification, install the EPI2ME application, select
`wf-metagenomics` at the version above, set the parameters in the table, and
point it at the pass reads. The equivalent Nextflow invocation is:

```bash
nextflow run epi2me-labs/wf-metagenomics -r v2.14.2 \
    --fastq fastq_pass --classifier kraken2 --database_set Standard-8 \
    --taxonomic_rank S --kraken2_confidence 0 --min_len 0 \
    --abundance_threshold 0 --bracken_threshold 10
```
