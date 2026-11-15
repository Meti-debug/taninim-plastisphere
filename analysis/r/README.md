# R analyses

Run from the repository root, in order. Each script sources `00_paths.R`,
which sets `DATA` to `data/` and `OUT` to `output/`.

```bash
Rscript analysis/r/03_beta_diversity_pcoa.R
```

| script | what it produces |
|---|---|
| `00_paths.R` | shared input and output paths; not run directly |
| `02_alpha_diversity.R` | Shannon diversity per sample type, Wilcoxon tests, boxplots |
| `03_beta_diversity_pcoa.R` | Bray-Curtis dissimilarity, PCoA, pairwise dissimilarity panels |
| `04_core_microbiome_heatmap.R` | the 15 most abundant species detected in all 12 samples of a sample type |
| `05_eskapee_heatmap.R` | ESKAPEE species abundance by site and sample type |
| `06_source_tracking_feast.R` | FEAST source tracking, four source-sink scenarios, and its stacked bars |
| `07_itol_annotations.R` | iTOL annotation files: branch colours by site, inner ring by sample type, outer ring by family |
| `08_dispersion_turnover.R` | betadisper dispersion tests and the within-site permutation test of river-sea turnover |
| `09_feast_sensitivity.R` | FEAST sensitivity for the S1 water sink: one or both plastic sources removed, a marine endmember added, and the full river pool, five seeds each |
| `10_river_signature.R` | non-marine taxa in river-mouth water: number, share of reads and the most abundant (Supplementary-6) |
| `11_feast_genus.R` | FEAST for the primary scenarios and the S1 water pools on genus-level counts (Supplementary-6) |
| `12_metric_robustness.R` | river or sea x sample type PERMANOVA and turnover test under six distance measures (Supplementary-1) |
| `13_marine_site_comparison.R` | river signature at the river mouth versus the reference beach, for plastic and water (Results 3.2) |

Inputs are described in the top-level `README.md`. `02` uses `sp_counts`, `sp_counts_clean` and
`meta_clean` built by `03`, so run `03` first in the same R session:

```r
source("analysis/r/03_beta_diversity_pcoa.R")
source("analysis/r/02_alpha_diversity.R")
``` `02`, `03` and `08` also need
`metadata.tsv`: one row per library, with site and sample type.

## Provenance

Scripts `02`–`07` are taken from the R session code used for the paper, with
absolute paths replaced by `DATA` and `OUT` and nothing else changed. `08` to
`13` were written afterwards, for the tests added in revision. `08` builds the Bray-Curtis
matrix exactly as `03` does. `09` needs the source-tracking FEAST from
github.com/cozygene/FEAST, not the Bioconductor package of the same name.

## Related code outside this folder

- The PERMANOVA reported in the manuscript was computed by
  `analysis/01_diversity_permanova.py`, with 9,999 permutations and sample type
  permuted within site. The `adonis2` calls in `02` and `03` use the package
  defaults (999 unrestricted permutations) and are not the reported test.
- The blaCTX-M group 1 amplicon phylogeny is in
  `analysis/phylogeny/ctxm/`, which has its own README: `run_all.sh`,
  `01_build_panel_and_tree.sh` for the allele panel, alignment and IQ-TREE run,
  `02_plot_complete_tree.R` and `03_plot_figure5c.R` for the figures. The alignment used for the tree, the
  treefile and the allele collapse table are in `output/ctxm/`.
