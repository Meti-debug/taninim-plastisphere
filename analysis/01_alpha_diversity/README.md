# 01 Alpha diversity

Shannon diversity of each library, and plastic-bound versus planktonic at each
site and across all sites.

Generates:
- Figure 3A-E (values and the significance marks)

Input:
- data/counts.csv

Output:
- analysis/01_alpha_diversity/output/alpha_diversity.json
- analysis/01_alpha_diversity/output/alpha_diversity.tsv

Run from the repository root:

```bash
python analysis/01_alpha_diversity/script.py
```

Method: Shannon index (natural log) on species relative abundances. Welch's
t-test per site (n = 3 per sample type) and across all sites (n = 12 per type);
Mann-Whitney U is reported alongside. See the header of `script.py`.
