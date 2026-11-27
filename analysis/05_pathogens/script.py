"""Pathogen abundance tests, plastic-associated vs planktonic (Figure 5A-B).

    python analysis/05_pathogens/script.py

INPUT
    data/counts.csv     Kraken2/Bracken species table (see README.md, Input tables)

OUTPUT (in analysis/05_pathogens/output/)
    fig5a_species_tests.tsv   per species: mean relative abundance (%)
                              by sample type, Mann-Whitney U p, BH q
    fig5b_eskapee_tests.tsv   total ESKAPEE abundance (CPM) by sample
                              type: Welch's t-test per site, Wilcoxon
                              rank-sum across all sites

METHOD
    Relative abundance and CPM are per library, against all classified reads.
    Fig. 5A  two-sided Mann-Whitney U, 12 plastic vs 12 water libraries per
             species; Benjamini-Hochberg correction across the species tested.
    Fig. 5B  ESKAPEE = E. faecium, S. aureus, K. pneumoniae, A. baumannii,
             P. aeruginosa, Enterobacter cloacae and E. coli, summed per library.
             Welch's t-test at each site (n = 3 per type) and Wilcoxon rank-sum
             (two-sided Mann-Whitney U) across all sites (n = 12 per type).
"""
import os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind

OUT = "analysis/05_pathogens/output"      # run from the repository root

SITES = ["R1", "R2", "S1", "Sc"]

# Species in Fig. 5A, in figure order.
FIG5A = ["Legionella clemsonensis", "Legionella pneumophila",
         "Legionella spiritensis", "Acinetobacter haemolyticus",
         "Aeromonas caviae", "Aeromonas dhakensis", "Aeromonas hydrophila",
         "Aeromonas media", "Aeromonas salmonicida", "Aeromonas veronii",
         "Klebsiella oxytoca", "Acinetobacter baumannii",
         "Klebsiella pneumoniae", "Pseudomonas aeruginosa",
         "Staphylococcus aureus", "Escherichia coli"]
ESKAPEE_SPECIES = ["Enterococcus faecium", "Staphylococcus aureus",
                   "Klebsiella pneumoniae", "Acinetobacter baumannii",
                   "Pseudomonas aeruginosa", "Enterobacter cloacae",
                   "Escherichia coli"]

# --- table, with sequencing labels converted to manuscript labels ------------
k = pd.read_csv("data/counts.csv")
libs = [c for c in k.columns if c[0] in "RS" and len(c) == 4]
# Sequencing labels (R2 = R1, R3 = R2, S1P4 = S1P1) are converted to manuscript
# labels; a table that already has manuscript labels is used as it is.
ren = dict((c, c) for c in libs)
if any(c.startswith("R3") for c in libs):
    for c in libs:
        n_ = "R1" + c[2:] if c.startswith("R2") else "R2" + c[2:] if c.startswith("R3") else c
        ren[c] = "S1P1" if n_ == "S1P4" else n_
X = k[libs].fillna(0).astype(float)
X.columns = [ren[c] for c in libs]
X = X[sorted(X.columns)]
assert X.shape[1] == 24, f"expected 24 libraries, found {X.shape[1]}"
species = k["species"].astype(str).str.strip()
total = X.sum(axis=0)
site = np.array([s[:2] for s in X.columns])
typ = np.array([s[2] for s in X.columns])
P, W = typ == "P", typ == "W"


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1)
    return out


# --- Fig. 5A ------------------------------------------------------------------
rows = []
for name in FIG5A:
    if not (species == name).any():
        raise SystemExit(f"{name} not found in counts.csv")
    ra = X[species == name].sum(axis=0) / total * 100
    a, b = ra[P].to_numpy(), ra[W].to_numpy()
    rows.append(dict(species=name, mean_plastic_pct=a.mean(), mean_water_pct=b.mean(),
                     mwu_p=mannwhitneyu(a, b, alternative="two-sided").pvalue))
fa = pd.DataFrame(rows)
fa["bh_q"] = bh(fa["mwu_p"])

# --- Fig. 5B ------------------------------------------------------------------
esk = species.isin(ESKAPEE_SPECIES)
cpm = (X[esk].sum(axis=0) / total * 1e6).to_numpy()
rows = []
for s in SITES:
    m = site == s
    a, b = cpm[m & P], cpm[m & W]
    rows.append(dict(comparison=s, test="Welch t", n=f"{len(a)} vs {len(b)}",
                     mean_plastic_cpm=a.mean(), mean_water_cpm=b.mean(),
                     p=ttest_ind(a, b, equal_var=False).pvalue))
a, b = cpm[P], cpm[W]
rows.append(dict(comparison="all sites", test="Wilcoxon rank-sum", n=f"{len(a)} vs {len(b)}",
                 mean_plastic_cpm=a.mean(), mean_water_cpm=b.mean(),
                 p=mannwhitneyu(a, b, alternative="two-sided").pvalue))
fb = pd.DataFrame(rows)

os.makedirs(OUT, exist_ok=True)
fa.to_csv(f"{OUT}/fig5a_species_tests.tsv", sep="\t", index=False)
fb.to_csv(f"{OUT}/fig5b_eskapee_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 140)
print(fa.round(4).to_string(index=False))
print()
print(fb.round(4).to_string(index=False))
