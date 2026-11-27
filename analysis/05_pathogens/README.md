# 05 Pathogens and ESKAPEE

Abundance of waterborne, opportunistic and ESKAPEE species on plastic versus in
water.

Generates:
- Figure 5A (significance marks per species)
- Figure 5B (ESKAPEE abundance and its tests)

Input:
- data/counts.csv

Output:
- analysis/05_pathogens/output/fig5a_species_tests.tsv
- analysis/05_pathogens/output/fig5b_eskapee_tests.tsv

Run from the repository root:

```bash
python analysis/05_pathogens/script.py
```

Method: Figure 5A, two-sided Mann-Whitney U, 12 plastic vs 12 water libraries
per species, Benjamini-Hochberg across the species shown. Figure 5B, ESKAPEE
counts per million summed per library; Welch's t-test per site (n = 3) and
Wilcoxon rank-sum across all sites (n = 12).
