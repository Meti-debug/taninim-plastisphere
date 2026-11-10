"""
Core taxa per habitat and sample type, and core retention from river to sea
(Fig. 4D).

Usage: python analysis/02_core_taxa.py

Definition (Fig. 4 caption): a species is core in a habitat when its relative
abundance is > 0.1% in at least 5 of the 6 samples of that habitat
(river = R1 + R2, marine = S1 + Sc; 3 replicates per site), computed
separately for plastic-bound (P) and planktonic (W) samples. Reads without a
species assignment ("Unknown") count towards the library total but are not
reported as a taxon.

Retention = share of river-core species that are also marine-core.
Test: is retention higher on plastic than in water? Observed statistic =
retention(plastic) - retention(water). Null model: sample type has no effect,
so within each site the six libraries (3 P + 3 W) are shuffled and three are
labelled plastic, three water; the statistic is recomputed 9,999 times.
One-sided p = (shuffles with statistic >= observed + 1) / (9,999 + 1).
"""
import numpy as np
import pandas as pd

THRESH, MIN_SAMPLES, N_PERM, SEED = 0.001, 5, 9999, 1
SITES = ["R1", "R2", "S1", "Sc"]
HAB = {"river": ["R1", "R2"], "marine": ["S1", "Sc"]}

k = pd.read_csv("data/counts.csv")
libs = [c for c in k.columns if c[0] in "RS" and len(c) == 4]
# Kraken2/Bracken table uses the older station labels (same as analysis/01):
# R2 -> R1, R3 -> R2, S1P4 -> S1P1.
ren = {}
for c in libs:
    n = "R1" + c[2:] if c.startswith("R2") else "R2" + c[2:] if c.startswith("R3") else c
    ren[c] = "S1P1" if n == "S1P4" else n
X = k[libs].fillna(0).rename(columns=ren)
X.index = k["species"].fillna("Unknown")
X = X.groupby(level=0).sum()                 # one row per species name
rel = X / X.sum(axis=0)                      # relative abundance per library
rel = rel.drop(index="Unknown", errors="ignore")
above = rel > THRESH                         # species x library, True/False
names = rel.index.to_numpy()


def core(cols):
    assert len(cols) == 6, cols
    return above[cols].sum(axis=1).to_numpy() >= MIN_SAMPLES


def summary(labels):
    """labels: dict site -> (3 plastic columns, 3 water columns)."""
    out = {}
    for t in (0, 1):
        r = core([c for s in HAB["river"] for c in labels[s][t]])
        m = core([c for s in HAB["marine"] for c in labels[s][t]])
        out[t] = (r, m)
    return out


def retention(r, m):
    return (r & m).sum() / r.sum()


obs_lab = {s: ([f"{s}P{i}" for i in (1, 2, 3)], [f"{s}W{i}" for i in (1, 2, 3)])
           for s in SITES}
obs = summary(obs_lab)
stat_obs = retention(*obs[0]) - retention(*obs[1])

for t, name in [(0, "Plastic-bound"), (1, "Planktonic")]:
    r, m = obs[t]
    print(f"{name}: river-core {r.sum()}, marine-core {m.sum()}, "
          f"core in both {(r & m).sum()} "
          f"(retention {100 * retention(r, m):.1f}% of river-core)")
    for s in sorted(names[r & m]):
        print(f"    {s}")

rng = np.random.default_rng(SEED)
hits = 0
for _ in range(N_PERM):
    lab = {}
    for s in SITES:
        six = rng.permutation(obs_lab[s][0] + obs_lab[s][1])
        lab[s] = (list(six[:3]), list(six[3:]))
    o = summary(lab)
    if retention(*o[0]) - retention(*o[1]) >= stat_obs - 1e-12:
        hits += 1
p = (hits + 1) / (N_PERM + 1)
print(f"Retention difference = {100 * stat_obs:.1f} percentage points; "
      f"one-sided within-site permutation p = {p:.4f} ({N_PERM} shuffles)")
