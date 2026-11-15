# The Taninim River-Sea Continuum: Plastisphere’s Persistence, Connectivity, and Potential Pathogen Dispersal into the Eastern Mediterranean Sea

Data analysis

Analysis code for the study.

Twenty-four libraries: 12 plastic-bound and 12 planktonic, 3 replicates at
each of 4 sites (R1 Ada Stream confluence, R2 estuary, S1 river mouth,
Sc Caesarea beach). Full-length 16S rRNA, Oxford Nanopore PromethION 2 Solo.

Two independent routes through the data, kept separate throughout:

- **Species-level profiles** from Kraken2 and Bracken, run through the EPI2ME
  `wf-metagenomics` workflow. Everything that needs a taxonomic name comes
  from here: diversity, ordination, PERMANOVA, core taxa, pathogen abundances,
  source tracking and the reference tree.
- **97% OTUs**, reference-free, built from cluster consensus sequences. Used
  only for OTU richness and for the shared and unique comparisons between
  sites and sample types. No OTU carries a taxonomic assignment and no Kraken2
  species is treated as an OTU; the two are never combined within a single
  statistic.

## Layout

```
preprocessing/   reads -> abundance tables (see its own README)
analysis/        Python statistics; writes intermediates to output/
figures/         Python plot scripts, and the CTX-M figure files
r/               R analyses and figures, one script per analysis
phylogeny/       reference 16S tree, and blaCTX-M group 1 amplicon tree in
                 phylogeny/ctxm/ (each has its own README)
data/            the input tables (not tracked; see data/README.md)
output/          everything generated; only output/ctxm/ is tracked (the
                 deposited CTX-M alignment, tree and IQ-TREE report)
```

Scripts use paths relative to the repository root, so run them from there:

```bash
mkdir -p output
python analysis/01_diversity_permanova.py
python figures/fig3_diversity_permanova.py
```

## Requirements

```
python >= 3.11        R >= 4.4
numpy >= 2.0          ggplot2 >= 3.5
pandas >= 2.0         pheatmap >= 1.0.12
scipy >= 1.11         dplyr, tidyr, scales
matplotlib >= 3.8
biopython >= 1.83     (phylogeny only)
```

External tools, for `phylogeny/` only: MAFFT v7, trimAl, IQ-TREE 2 v2.3.5.
The OTU pipeline needs VSEARCH v2.32.0.

Figure scripts use the Liberation Sans font.

## Which script produces what

`analysis/` writes the intermediates that `figures/` reads, so analysis comes
first. Within each group the numeric prefix is the order.

### preprocessing

| script | what it does |
|---|---|
| `01_epi2me_kraken2.md` | record of the EPI2ME run: workflow version, database and every parameter |
| `02_epi2me_counts_to_table.py` | EPI2ME species counts -> `data/counts.csv` |
| `03_otu97_pipeline.sh` | quality filter, consensus clustering, chimera removal, mapping -> `otus.fasta`, `otu_table.tsv` |

### analysis

| script | what it produces | reads |
|---|---|---|
| `01_diversity_permanova.py` | Shannon, Bray–Curtis, PCoA, PERMANOVA -> `output/fig3.json` | `data/counts.csv` |
| `02_core_taxa.py` | core taxa per habitat and sample type, and core retention from river to sea with its permutation test | `data/counts.csv` |
| `03_four_site_persistence.py` | OTUs present at all four sites, plastic vs planktonic, with its permutation test | an `otu_table.tsv` |
| `04_discovery_depth_sensitivity.py` | side-by-side comparison of two OTU tables built at different discovery depths | two `otu_table.tsv` |

### figures

| script | what it draws | reads |
|---|---|---|
| `fig3_diversity_permanova.py` | diversity, ordination, Bray–Curtis structure | `01` |
| `fig4ab_venn_otu.py` | shared and unique OTUs, by site and by sample type | an `otu_table.tsv` |

### r

Run in order from the repository root; each sources `r/00_paths.R`.

| script | what it produces |
|---|---|
| `01_chromagar_colonies.R` | CHROMagar colony colours and total growth per river |
| `02_alpha_diversity.R` | Shannon diversity per sample type, Wilcoxon tests |
| `03_beta_diversity_pcoa.R` | Bray-Curtis, PCoA, pairwise dissimilarity panels |
| `04_core_microbiome_heatmap.R` | the 15 most abundant consistently detected species |
| `05_eskapee_heatmap.R` | ESKAPEE species abundance by site and sample type |
| `06_source_tracking_feast.R` | FEAST source tracking and its stacked bars |
| `07_itol_annotations.R` | iTOL annotation rings for the reference tree |

`figures_as_run.R` is the full session code, kept verbatim as the provenance
record. See `r/README.md`.

### phylogeny

`bash phylogeny/run_all.sh` builds the reference 16S tree from the 500 most
abundant species. See `phylogeny/README.md`.

### CTX-M group 1 tree (Figure 5C and the complete tree)

The tree is deposited in `output/ctxm/`, so the figures can be redrawn without
re-running IQ-TREE. From the repository root:

```bash
Rscript phylogeny/ctxm/02_plot_complete_tree.R   # figures/Figure_CTXM_sample_tree.pdf/.png
Rscript phylogeny/ctxm/03_plot_figure5c.R        # figures/Figure5C_CTXM.pdf/.png
```

| file | role |
|---|---|
| `phylogeny/ctxm/01_build_panel_and_tree.sh` | curate references, align (MAFFT), cut to CDS 28-400, collapse identical alleles, infer the tree (IQ-TREE) |
| `phylogeny/ctxm/02_plot_complete_tree.R` | complete tree, all tips, SH-aLRT/UFBoot at every resolved node |
| `phylogeny/ctxm/03_plot_figure5c.R` | Figure 5C, the tree pruned to the tips in `figure5c_tips.txt` |
| `phylogeny/ctxm/plot_helpers.R` | shared plotting functions |
| `output/ctxm/ctxm_group1.treefile` | the tree both figures are drawn from |
| `output/ctxm/ctxm_group1.iqtree`, `.log` | model selection, alignment statistics, support |
| `output/ctxm/tree_input.fasta` | the alignment the tree was inferred from |
| `output/ctxm/haplotype_collapse.tsv` | which alleles each collapsed reference tip stands for |

R packages: ggtree, treeio, ape, ggplot2, dplyr (see `phylogeny/ctxm/README.md`
for a conda command). Re-running `01` needs the raw consensus and reference
sequences, `export CTXM_DATA=...`, then `bash phylogeny/ctxm/run_all.sh`.

## Reproducing the OTU table

```bash
IN_DIR=path/to/fastq OUT_DIR=otu97_25k BUILD_N=25000 \
  bash preprocessing/03_otu97_pipeline.sh
```

`BUILD_N` is the number of quality-filtered reads per library pooled to
construct the OTU set. The reported analysis uses 25,000; the sensitivity run
uses 10,000. Every other parameter is identical between the two, so discovery
depth is the only variable. `04_discovery_depth_sensitivity.py` compares them.

## Data

Raw reads: NCBI SRA, BioProject PRJNA1458552. The input tables are described
in `data/README.md` and deposited with the manuscript as Supplementary Data.
