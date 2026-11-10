"""
Four-site persistent OTUs: plastic-bound vs planktonic.

Usage: python analysis/03_four_site_persistence.py OTU_TABLE

Presence rule (same as Fig. 4A/B): an OTU is present at a site, in one
sample type, when it has >= 8 reads in >= 2 of the 3 replicates.

Statistic: the share of OTUs that are present at all four sites
(R1, R2, S1, Sc), computed separately for plastic and for water;
observed difference = plastic share - water share.

Null model: sample type has no effect. Within each site the six
libraries (3 P + 3 W) are shuffled and three are labelled "plastic",
the other three "water". The statistic is recomputed for every shuffle.
One-sided p = (number of shuffles with difference >= observed + 1) / (N + 1).
"""
import sys
import numpy as np
import pandas as pd

MIN_READS, MINPREV, N_PERM, SEED = 8, 2, 9999, 1
SITES = ["R1", "R2", "S1", "Sc"]

tab = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
hit = (tab >= MIN_READS).to_numpy()          # OTU x library, True/False
cols = list(tab.columns)


def share_all_sites(groups):
    """groups: dict site -> list of 3 column indices. Returns (n_all4, n_total)."""
    pres = np.stack([hit[:, idx].sum(axis=1) >= MINPREV for idx in groups.values()],
                    axis=1)                    # OTU x site
    n_total = pres.any(axis=1).sum()
    n_all4 = pres.all(axis=1).sum()
    return n_all4, n_total


def stat(p_groups, w_groups):
    a_p, t_p = share_all_sites(p_groups)
    a_w, t_w = share_all_sites(w_groups)
    return a_p / t_p - a_w / t_w, (a_p, t_p, a_w, t_w)


obs_p = {s: [cols.index(f"{s}P{i}") for i in (1, 2, 3)] for s in SITES}
obs_w = {s: [cols.index(f"{s}W{i}") for i in (1, 2, 3)] for s in SITES}
obs, (a_p, t_p, a_w, t_w) = stat(obs_p, obs_w)

rng = np.random.default_rng(SEED)
null = np.empty(N_PERM)
for k in range(N_PERM):
    pg, wg = {}, {}
    for s in SITES:
        six = rng.permutation(obs_p[s] + obs_w[s])
        pg[s], wg[s] = list(six[:3]), list(six[3:])
    null[k] = stat(pg, wg)[0]

p = (np.sum(null >= obs) + 1) / (N_PERM + 1)
print(f"Plastic: {a_p} of {t_p} OTUs at all four sites ({100*a_p/t_p:.2f}%)")
print(f"Water:   {a_w} of {t_w} OTUs at all four sites ({100*a_w/t_w:.2f}%)")
print(f"Difference = {100*obs:.2f} percentage points; "
      f"one-sided permutation p = {p:.4f} ({N_PERM} within-site shuffles)")
