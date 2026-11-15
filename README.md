# The Taninim River-Sea Continuum: Plastisphere’s Persistence, Connectivity, and Potential Pathogen Dispersal into the Eastern Mediterranean Sea

Data analysis

Analysis code for the study.

- **Species-level profiles** from Kraken2 and Bracken, run through the EPI2ME
  `wf-metagenomics` workflow. Everything that needs a taxonomic name comes
  from here: diversity, ordination, PERMANOVA, core taxa, pathogen abundances,
  source tracking and the reference tree.
- **97% OTUs**, reference-free, built from cluster consensus sequences. Used
  only for OTU richness and for the shared and unique comparisons between
  sites and sample types. No OTU carries a taxonomic assignment and no Kraken2
  species is treated as an OTU; the two are never combined within a single
  statistic.


## Requirements

```
python >= 3.11        R >= 4.4
numpy >= 2.0          ggplot2 >= 3.5
pandas >= 2.0         pheatmap >= 1.0.12
scipy >= 1.11         dplyr, tidyr, scales
biopython >= 1.83     (phylogeny only)
```

External tools, for `analysis/phylogeny/` only: MAFFT v7, trimAl, IQ-TREE 2 v2.3.5.

## Which script produces what

Within each group the numeric prefix is the order.

### analysis

| script | what it produces | reads |
|---|---|---|
| `01_diversity_permanova.py` | Shannon, Bray–Curtis, PCoA, PERMANOVA -> `output/fig3.json` | `data/counts.csv` |
| `02_core_taxa.py` | core taxa per habitat and sample type, and core retention from river to sea with its permutation test | `data/counts.csv` |
| `03_four_site_persistence.py` | OTUs present at all four sites, plastic vs planktonic, with its permutation test | an `otu_table.tsv` |
| `04_discovery_depth_sensitivity.py` | side-by-side comparison of two OTU tables built at different discovery depths | two `otu_table.tsv` |
| `05_pathogen_tests.py` | Fig. 5A Mann–Whitney U with BH correction per species; Fig. 5B ESKAPEE CPM, Welch per site and Wilcoxon rank-sum overall | `data/counts.csv` |

### analysis/r

FEAST source tracking and the dispersion and robustness tests, which are run in
R. Each script runs on its own; see `analysis/r/README.md`.

| script | what it produces |
|---|---|
| `01_source_tracking_feast.R` | FEAST source tracking, four source-sink scenarios |
| `02_dispersion_turnover.R` | betadisper dispersion and river-sea turnover tests |
| `03_feast_sensitivity.R` | FEAST sensitivity for the S1 water sink |
| `04_river_signature.R` | non-marine taxa in river-mouth water |
| `05_feast_genus.R` | FEAST on genus-level counts |
| `06_metric_robustness.R` | PERMANOVA and turnover under six distance measures |
| `07_marine_site_comparison.R` | river signature at the river mouth versus the reference beach |

### analysis/phylogeny

`bash analysis/phylogeny/run_all.sh` builds the reference 16S tree from the 500 most
abundant species. See `analysis/phylogeny/README.md`.

### CTX-M group 1 tree (Figure 5C and the complete tree)

The tree is deposited in `output/ctxm/`, so the figures can be redrawn without
re-running IQ-TREE. From the repository root:

```bash
Rscript analysis/phylogeny/ctxm/02_plot_complete_tree.R   # output/figures/Figure_CTXM_sample_tree.pdf/.png
Rscript analysis/phylogeny/ctxm/03_plot_figure5c.R        # output/figures/Figure5C_CTXM.pdf/.png
```

Panel C as placed in the composite Figure 5 (1758 x 1044 px at 600 dpi, text
scaled to match panels A and B):

```bash
CTXM_STEM=Figure5C_CTXM_composite CTXM_W=2.93 CTXM_H=1.74 CTXM_TEXT=0.66 \
CTXM_TREE_FR=0.40 CTXM_ROOT_LEN=0.013 Rscript analysis/phylogeny/ctxm/03_plot_figure5c.R
```

| file | role |
|---|---|
| `analysis/phylogeny/ctxm/01_build_panel_and_tree.sh` | curate references, align (MAFFT), cut to CDS 28-400, collapse identical alleles, infer the tree (IQ-TREE) |
| `analysis/phylogeny/ctxm/02_plot_complete_tree.R` | complete tree, all tips, SH-aLRT/UFBoot at every resolved node |
| `analysis/phylogeny/ctxm/03_plot_figure5c.R` | Figure 5C, the tree pruned to the tips in `figure5c_tips.txt` |
| `analysis/phylogeny/ctxm/plot_helpers.R` | shared plotting functions |
| `output/ctxm/ctxm_group1.treefile` | the tree both figures are drawn from |
| `output/ctxm/ctxm_group1.iqtree`, `.log` | model selection, alignment statistics, support |
| `output/ctxm/tree_input.fasta` | the alignment the tree was inferred from |
| `output/ctxm/haplotype_collapse.tsv` | which alleles each collapsed reference tip stands for |

R packages: ggtree, treeio, ape, ggplot2, dplyr (see `analysis/phylogeny/ctxm/README.md`
for a conda command). Re-running `01` needs the raw consensus and reference
sequences, `export CTXM_DATA=...`, then `bash analysis/phylogeny/ctxm/run_all.sh`.

## Input tables

The scripts read two tables from a `data/` folder in the repository root.
Both are deposited with the manuscript as Supplementary Data; raw reads are in
the NCBI SRA under BioProject PRJNA1458552.

| file | contents |
|---|---|
| `data/counts.csv` | Kraken2/Bracken species counts from EPI2ME `wf-metagenomics`: `tax`, `genus`, `species`, then one column per library |
| `data/otu_table.tsv` | 97% OTU x library read counts, 24 libraries |

`counts.csv` carries the site labels used during sequencing, one position
upstream of the manuscript labels; the scripts convert them on read. The OTU
table uses manuscript labels.

| in `counts.csv` | in the manuscript |
|---|---|
| `R2P1`–`R2W3` | `R1P1`–`R1W3` |
| `R3P1`–`R3W3` | `R2P1`–`R2W3` |
| `S1P4` | `S1P1` |
| all others | unchanged |
