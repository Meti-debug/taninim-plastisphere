# R analyses

Run from the repository root, in order. Each script sources `00_paths.R`,
which sets `DATA` to `data/` and `OUT` to `output/`.

```bash
Rscript r/01_chromagar_colonies.R
```

| script | what it produces |
|---|---|
| `00_paths.R` | shared input and output paths; not run directly |
| `01_chromagar_colonies.R` | CHROMagar colony colours and total growth per river |
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

Inputs are described in `data/README.md`. `02`, `03` and `08` also need
`metadata.tsv`: one row per library, with site and sample type.

## Provenance

`figures_as_run.R` is the full R session code, deposited verbatim. Scripts
`01`–`07` are extracted from it, with absolute paths replaced by
`DATA` and `OUT` and nothing else changed; each script's header gives the line
range it came from. The session file also contains superseded versions of some
figures, a network analysis and Venn panels that are not in the manuscript, and
a Colab notebook block that built the reference 16S tree. It is kept as the
record of what was run.

`08` to `13` were written after the
session file, for the tests added in revision. `08` builds the Bray-Curtis
matrix exactly as `03` does. `09` needs the source-tracking FEAST from
github.com/cozygene/FEAST, not the Bioconductor package of the same name.

## Not in this folder

- The PERMANOVA reported in the manuscript was computed by
  `analysis/01_diversity_permanova.py`, with 9,999 permutations and sample type
  permuted within site. The `adonis2` calls in `figures_as_run.R` use the
  package defaults of 999 unrestricted permutations and are not the reported
  test.
- The blaCTX-M group 1 amplicon phylogeny is not in this folder. It is in
  `phylogeny/ctxm/`, which has its own README: `run_all.sh`,
  `01_build_panel_and_tree.sh` for the allele panel, alignment and IQ-TREE run,
  and `02` to `06` for the figure panels. The alignment used for the tree, the
  treefile and the allele collapse table are in `output/ctxm/`.
