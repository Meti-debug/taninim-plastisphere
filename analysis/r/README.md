# R analyses

Source tracking (FEAST) and the dispersion and robustness tests. Each script
runs on its own from the repository root and sources `00_paths.R`, which sets
`DATA` to `data/` and `OUT` to `output/`.

```bash
Rscript analysis/r/01_source_tracking_feast.R
```

| script | what it produces |
|---|---|
| `00_paths.R` | shared input and output paths; not run directly |
| `01_source_tracking_feast.R` | FEAST source tracking for the four source-sink scenarios (Fig. 4E) |
| `02_dispersion_turnover.R` | betadisper dispersion tests and the within-site permutation test of river-sea turnover (Table S1, Results 3.2) |
| `03_feast_sensitivity.R` | FEAST for the S1 water sink with one or both plastic sources removed, a marine endmember added, and the full river pool, five seeds each (Table S6, Fig. S7) |
| `04_river_signature.R` | non-marine taxa in river-mouth water: number, share of reads and the most abundant (Table S6) |
| `05_feast_genus.R` | FEAST for the primary scenarios and the S1 water pools on genus-level counts (Table S6) |
| `06_metric_robustness.R` | river or sea x sample type PERMANOVA and turnover test under six distance measures (Table S1) |
| `07_marine_site_comparison.R` | river signature at the river mouth versus the reference beach, for plastic and water (Results 3.2) |

Inputs are the tables described in the top-level `README.md`. `01` and `02`
also need `data/metadata.tsv`: one row per library, with site and sample type.

FEAST is the source-tracking package from github.com/cozygene/FEAST (v0.1.0),
not the Bioconductor package of the same name.

The PERMANOVA reported in the paper is computed by
`analysis/01_diversity_permanova.py` (9,999 permutations, sample type permuted
within site). The `adonis2` calls in `02` and `06` are the river or sea x
sample type checks reported in the supplementary tables.
