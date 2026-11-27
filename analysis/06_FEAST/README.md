# 06 FEAST source tracking

Source tracking from river to sea, its sensitivity analyses, and the river
signature in river-mouth water.

Generates:
- Figure 4E (source contributions for four sinks)
- Table S6 and Figure S7 (FEAST sensitivity, genus level, river signature)
- Results 3.2 comparison of the river mouth with the reference beach

Input:
- data/counts.csv

Output (analysis/06_FEAST/output/):
- feast_primary/feast_primary_summary.csv (mean and SD per sink and source), feast_primary_all_runs.csv, Figure4E_FEAST.pdf/.png
- feast_sensitivity/feast_sensitivity_summary.csv
- feast_genus/feast_genus_summary.csv
- `river_signature.R` and `marine_site_comparison.R` print to the console

Run from the repository root:

```bash
Rscript analysis/06_FEAST/primary_analysis.R
Rscript analysis/06_FEAST/sensitivity_analysis.R
Rscript analysis/06_FEAST/genus_level.R
Rscript analysis/06_FEAST/river_signature.R
Rscript analysis/06_FEAST/marine_site_comparison.R
```

| script | what it does |
|---|---|
| `primary_analysis.R` | FEAST for the four primary source-sink scenarios, rarefied to 91,000 reads, seeds 1-5, and the stacked bars of Figure 4E |
| `sensitivity_analysis.R` | the S1 water sink with eight different source pools, rarefied to 91,000 reads, seeds 1-5 |
| `genus_level.R` | the primary scenarios and the S1 water pools on genus-level counts |
| `river_signature.R` | non-marine taxa in river-mouth water: number, share of reads, the most abundant |
| `marine_site_comparison.R` | river signature at the river mouth versus the reference beach (exact permutation, 20 splits) |

All FEAST scripts filter species the same way (>= 50 reads in total and
present in >= 2 replicates at one site), rarefy every library to 91,000 reads
with seed 42, run each scenario with seeds 1 to 5, and report the mean and SD
over 3 sink replicates x 5 runs.

FEAST is the package from github.com/cozygene/FEAST (v0.1.0), not the
Bioconductor package of the same name. FEAST changes the working directory
while it runs, so output paths are absolute (`analysis/paths.R`).
