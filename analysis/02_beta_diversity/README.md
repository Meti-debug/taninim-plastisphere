# 02 Beta diversity

Bray-Curtis dissimilarity, PCoA and PERMANOVA, and the tests of dispersion and
river-sea turnover.

Generates:
- Figure 3F (PCoA, PERMANOVA R² for site and sample type)
- Figure 3G (within- and between-group dissimilarity, compact letters)
- Figure 3H (pairwise site dissimilarity by sample type)
- Table S1 (dispersion, turnover, PERMANOVA under six distance measures)

Input:
- data/counts.csv
- data/metadata.tsv (R scripts only)

Output:
- analysis/02_beta_diversity/output/beta_diversity.json
- analysis/02_beta_diversity/output/metric_robustness.csv
- `dispersion_turnover.R` prints its tables to the console

Run from the repository root:

```bash
python  analysis/02_beta_diversity/script.py
Rscript analysis/02_beta_diversity/dispersion_turnover.R
Rscript analysis/02_beta_diversity/metric_robustness.R
```

| script | what it does |
|---|---|
| `script.py` | Bray-Curtis on read counts, PCoA, PERMANOVA (9,999 permutations, sample type permuted within site), panel G letters, panel H distances. Add `relative` as an argument for Bray-Curtis on proportions. |
| `dispersion_turnover.R` | betadisper dispersion tests, within-site permutation test of river-sea turnover, river or sea x sample type PERMANOVA |
| `metric_robustness.R` | the PERMANOVA and turnover test repeated under six distance measures |

The PERMANOVA reported in the paper is the one in `script.py`. The `adonis2`
calls in the R scripts are the river or sea x sample type checks in Table S1.
