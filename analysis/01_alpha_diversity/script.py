"""Alpha diversity: Shannon index per library and per site (Figure 3A-E).

    python analysis/01_alpha_diversity/script.py

INPUT
    data/counts.csv     Kraken2/Bracken species table (see data/README.md).

OUTPUT (in analysis/01_alpha_diversity/output/)
    alpha_diversity.json   per-library Shannon index, and for each site and
                           for all sites together the plastic and water
                           values with Welch's t-test and Mann-Whitney U p.
    alpha_diversity.tsv    the per-library Shannon index as a table.

METHOD
    Shannon     natural log, on species relative abundances, zeros dropped.
    Tests       Welch's t-test, plastic vs water, at each site (n = 3 per type,
                Fig. 3A-D) and across all sites (n = 12 per type, Fig. 3E).
                Mann-Whitney U is reported alongside.  Within one site there
                are only 20 ways to split six libraries into three plastic and
                three water, so a permutation test cannot give a two-sided
                p below 0.10; this is why the per-site test is Welch's t-test.
"""
import os
import json
import numpy as np, pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu

OUT = "analysis/01_alpha_diversity/output"      # run from the repository root
SITES = ["R1", "R2", "S1", "Sc"]

# The Kraken2/Bracken table carries the station labels used during
# sequencing: R2 -> R1, R3 -> R2, and S1P4 -> S1P1.
k = pd.read_csv("data/counts.csv")
ks = [c for c in k.columns if c[0] in "RS" and len(c) == 4]
# Sequencing labels (R2 = R1, R3 = R2, S1P4 = S1P1)
ren = dict((c, c) for c in ks)
if any(c.startswith("R3") for c in ks):
    for c in ks:
        n_ = "R1" + c[2:] if c.startswith("R2") else "R2" + c[2:] if c.startswith("R3") else c
        ren[c] = "S1P1" if n_ == "S1P4" else n_
K = k[ks].fillna(0).copy(); K.columns = [ren[c] for c in ks]
sm = sorted(K.columns)
assert len(sm) == 24, f"expected 24 libraries, found {len(sm)}"
rel = K[sm].to_numpy(float).T
rel = rel / rel.sum(1, keepdims=True)
site = np.array([s[:2] for s in sm]); typ = np.array([s[2] for s in sm])
n = len(sm)

H = np.array([-(p[p > 0] * np.log(p[p > 0])).sum() for p in rel])

# ---------- alpha per site -------------------------------------------------
alpha = {}
for s_ in SITES + ["ALL"]:
    m = np.ones(n, bool) if s_ == "ALL" else (site == s_)
    a, b = H[m & (typ == "P")], H[m & (typ == "W")]
    alpha[s_] = dict(P=a.tolist(), W=b.tolist(),
                     t_p=float(ttest_ind(a, b, equal_var=False).pvalue),
                     mwu_p=float(mannwhitneyu(a, b, alternative="two-sided").pvalue))
    print(f"{s_:4s} plastic {a.mean():.2f} water {b.mean():.2f}  Welch p={alpha[s_]['t_p']:.2e}")

os.makedirs(OUT, exist_ok=True)
json.dump(dict(sample=sm, site=site.tolist(), type=typ.tolist(),
               shannon=H.tolist(), alpha=alpha),
          open(f"{OUT}/alpha_diversity.json", "w"), indent=1)
pd.DataFrame(dict(sample=sm, site=site, type=typ, shannon=H)).to_csv(
    f"{OUT}/alpha_diversity.tsv", sep="\t", index=False, float_format="%.4f")
print("saved")
