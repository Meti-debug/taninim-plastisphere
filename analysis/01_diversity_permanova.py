"""Alpha diversity, beta diversity, ordination and PERMANOVA (Figure 3).

WHAT IT DOES
    The whole community-structure analysis, from the species table to the
    numbers behind every panel of Figure 3.  Nothing here calls vegan or
    scikit-bio; the four methods are implemented directly so that the exact
    model being fitted is visible, and so that the permutation scheme can be
    restricted the way this design requires.

    python analysis/01_diversity_permanova.py
    python analysis/01_diversity_permanova.py raw   # Bray-Curtis on raw counts

INPUT
    data/counts.csv             Kraken2/Bracken species table, written by
                                preprocessing/02_epi2me_counts_to_table.py.

OUTPUT
    output/fig3.json            everything the figure script needs: per-sample
                                Shannon, PCoA coordinates and explained
                                variance, PERMANOVA table, and the
                                within/between-group Bray-Curtis distributions
                                for panels G and H.

METHOD
    Shannon          natural log, on species relative abundances, zeros dropped.
    Bray-Curtis      sum|p_i - p_j| / sum(p_i + p_j) on relative abundances, so
                     the 95,525-read library is not penalised for its depth.
    PCoA             double-centred Gower matrix G = -0.5 * C D^2 C, symmetric
                     eigendecomposition, coordinates scaled by sqrt(eigenvalue);
                     negative eigenvalues are discarded, not corrected.
    PERMANOVA        sequential (adonis2-style) sums of squares from traces of
                     hat matrices, trace(H G H), with 9,999 permutations.
                     Two models are fitted, matching the two reported in the
                     manuscript.  The first is station then sample type, and
                     gives the station and sample-type terms.  The second is
                     river-or-sea, then sample type, then their interaction,
                     and gives the interaction term only.
                     For the sample-type term, permutations are RESTRICTED to
                     within station: the three plastic and three water
                     replicates are exchangeable with each other but not
                     across stations, and an unrestricted permutation would
                     test a null nobody is proposing.
                     The interaction term is currently permuted WITHOUT that
                     restriction, and the permuted labels are those of the
                     combined factor rather than of sample type.  The
                     restricted scheme is the one the design argues for, and
                     this term is to be moved onto it.

READ THIS BEFORE QUOTING A p-VALUE
    With four stations x two types x three replicates the permutable space is
    small.  A restricted permutation test on sample type has 20^4 = 160,000
    distinct arrangements, which is plenty; a per-station test has only 20, so
    the smallest two-sided p obtainable at one station is 0.10.  That floor is
    a property of the design, not of the software.
"""
import numpy as np, pandas as pd, itertools, json
from scipy.stats import ttest_ind, mannwhitneyu

rng = np.random.default_rng(0)
SITES = ["R1", "R2", "S1", "Sc"]

import sys
# The Kraken2/Bracken table carries the station labels used during
# sequencing: R2 -> R1, R3 -> R2, and S1P4 -> S1P1.  The mapping is recorded
# in data/README.md.
k = pd.read_csv("data/counts.csv")
ks = [c for c in k.columns if c[0] in "RS" and len(c) == 4]
ren = {}
for c in ks:
    n_ = "R1" + c[2:] if c.startswith("R2") else "R2" + c[2:] if c.startswith("R3") else c
    ren[c] = "S1P1" if n_ == "S1P4" else n_
K = k[ks].fillna(0).copy(); K.columns = [ren[c] for c in ks]
sm = sorted(K.columns)
assert len(sm) == 24, f"expected 24 libraries, found {len(sm)}"
rel = K[sm].to_numpy(float).T
SOURCE = "Kraken2/Bracken"
print(f"classifier: {SOURCE}")
# Bray-Curtis input. "relative" (default): proportions, so library depth does
# not enter the distance. "raw": read counts as they are, which is what
# vegan::vegdist(t(counts), "bray") computes and what the original Figure 3
# used; with unequal library depth it inflates the distance between shallow
# and deep libraries.
BC_INPUT = "raw" if "raw" in sys.argv[1:] else "relative"
counts_ = rel.copy()
rel = rel / rel.sum(1, keepdims=True)
X_BC = counts_ if BC_INPUT == "raw" else rel
print(f"Bray-Curtis on {BC_INPUT} abundances")
site = np.array([s[:2] for s in sm]); typ = np.array([s[2] for s in sm])
hab = np.where(np.isin(site, ["R1", "R2"]), "river", "marine")
n = len(sm)

H = np.array([-(p[p > 0] * np.log(p[p > 0])).sum() for p in rel])

BC = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        BC[i, j] = BC[j, i] = (np.abs(X_BC[i] - X_BC[j]).sum()
                               / (X_BC[i] + X_BC[j]).sum())

# ---------- PCoA -----------------------------------------------------------
A = -0.5 * BC ** 2
C = np.eye(n) - np.ones((n, n)) / n
G = C @ A @ C
w, V = np.linalg.eigh(G); o = np.argsort(w)[::-1]; w, V = w[o], V[:, o]
pos = w > 1e-9
CO = V[:, pos] * np.sqrt(w[pos]); EXP = w[pos] / w[pos].sum()

# ---------- PERMANOVA (adonis2-style, sequential) --------------------------
def design(*factors):
    cols = [np.ones((n, 1))]
    for f in factors:
        lv = np.unique(f)[1:]
        cols.append(np.column_stack([(f == L).astype(float) for L in lv]))
    return np.hstack(cols)

def hat(*factors):
    X = design(*factors)
    return X @ np.linalg.pinv(X.T @ X) @ X.T

def term_dfs(factor_list):
    """Degrees of freedom each factor adds, in order, as the increase in the
    rank of the design matrix. An interaction coded as a combined factor adds
    only the columns not already spanned by its main effects."""
    out, prev = [], 1
    for k in range(1, len(factor_list) + 1):
        r = np.linalg.matrix_rank(design(*factor_list[:k]))
        out.append(r - prev); prev = r
    return out

def ss(Hm):
    return float(np.trace(Hm @ G @ Hm))

SST = float(np.trace(G))

def seq_terms(factor_list):
    """Sequential SS for each factor in order."""
    out, prev = [], 0.0
    for k in range(1, len(factor_list) + 1):
        cur = ss(hat(*factor_list[:k]))
        out.append(cur - prev); prev = cur
    return out, prev

def permanova_term(factor_list, which, restrict=None, nperm=9999):
    """p for the `which`-th factor, permuting that factor's labels."""
    terms, ssmodel = seq_terms(factor_list)
    ssres = SST - ssmodel
    dfs = term_dfs(factor_list)
    dfres = n - 1 - sum(dfs)
    Fobs = (terms[which] / dfs[which]) / (ssres / dfres)
    hits = 0
    for _ in range(nperm):
        fl = [f.copy() for f in factor_list]
        tgt = fl[which]
        if restrict is None:
            fl[which] = rng.permutation(tgt)
        else:
            g = tgt.copy()
            for s_ in np.unique(restrict):
                k = np.where(restrict == s_)[0]
                g[k] = rng.permutation(tgt[k])
            fl[which] = g
        t2, m2 = seq_terms(fl)
        F2 = (t2[which] / dfs[which]) / ((SST - m2) / dfres)
        if F2 >= Fobs - 1e-12:
            hits += 1
    return dict(R2=terms[which] / SST, F=float(Fobs), p=(hits + 1) / (nperm + 1))

pf_site = permanova_term([site, typ], 0, restrict=None, nperm=9999)
pf_type = permanova_term([site, typ], 1, restrict=site, nperm=9999)
print("panel F  site :", {k: round(v, 4) for k, v in pf_site.items()})
print("panel F  type :", {k: round(v, 4) for k, v in pf_type.items()})

ph_int = permanova_term([hab, typ, np.char.add(np.char.add(hab, "_"), typ)], 2,
                        restrict=None, nperm=9999)
print("panel H  habitat x type :", {k: round(v, 4) for k, v in ph_int.items()})

# ---------- alpha per site -------------------------------------------------
alpha = {}
for s_ in SITES + ["ALL"]:
    m = np.ones(n, bool) if s_ == "ALL" else (site == s_)
    a, b = H[m & (typ == "P")], H[m & (typ == "W")]
    alpha[s_] = dict(P=a.tolist(), W=b.tolist(),
                     t_p=float(ttest_ind(a, b, equal_var=False).pvalue),
                     mwu_p=float(mannwhitneyu(a, b, alternative="two-sided").pvalue))
    print(f"{s_:4s} plastic {a.mean():.2f} water {b.mean():.2f}  Welch p={alpha[s_]['t_p']:.2e}")

# ---------- panel G --------------------------------------------------------
GK = ["P-P\nper site", "P-P\ndiff site", "W-W\nper site", "W-W\ndiff site",
      "P-W\nper site", "P-W\ndiff site"]
gg = {k: [] for k in GK}
for i, j in itertools.combinations(range(n), 2):
    same = site[i] == site[j]
    base = "P-P" if typ[i] == typ[j] == "P" else "W-W" if typ[i] == typ[j] == "W" else "P-W"
    gg[f"{base}\n{'per' if same else 'diff'} site"].append(BC[i, j])

# ---------- panel G letters ------------------------------------------------
# Same test as the original R analysis (multcompView): the statistic is the
# MEDIAN Bray-Curtis dissimilarity of each of the six comparison classes, and
# the null distribution comes from permuting the SAMPLE labels (sample type and
# site together) 9,999 times while the distance matrix stays fixed. Permuting
# samples, not distances, respects the fact that every sample takes part in
# many pairs. Two-sided p for each pair of classes, Holm-adjusted, then a
# compact letter display (insert-and-absorb, Piepho 2004): classes sharing a
# letter do not differ at p_holm < 0.05.
N_PERM_G = 9999
iu, ju = np.triu_indices(n, 1)
D_pairs = BC[iu, ju]

def pair_class(t_, s_):
    same = s_[iu] == s_[ju]
    base = np.where((t_[iu] == "P") & (t_[ju] == "P"), 0,
                    np.where((t_[iu] == "W") & (t_[ju] == "W"), 2, 4))
    return base + (~same).astype(int)          # GK order: PP same, PP diff, ...

def class_medians(t_, s_):
    c = pair_class(t_, s_)
    return np.array([np.median(D_pairs[c == k]) if np.any(c == k) else np.nan
                     for k in range(6)])

obs_med = class_medians(typ, site)
prng = np.random.default_rng(1)
null_med = np.empty((N_PERM_G, 6))
for b in range(N_PERM_G):
    idx = prng.permutation(n)
    null_med[b] = class_medians(typ[idx], site[idx])

def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    adj = np.maximum.accumulate((m - np.arange(m)) * p[o])
    out = np.empty(m); out[o] = np.minimum(adj, 1.0)
    return out

g_tests = []
for a_, b_ in itertools.combinations(range(6), 2):
    obs = obs_med[b_] - obs_med[a_]
    nul = null_med[:, b_] - null_med[:, a_]
    nul = nul[~np.isnan(nul)]
    p_ = (np.sum(np.abs(nul) >= abs(obs)) + 1) / (len(nul) + 1)
    g_tests.append(dict(a=GK[a_], b=GK[b_], test="label permutation, median",
                        p=float(p_)))
for t_, q_ in zip(g_tests, holm([t["p"] for t in g_tests])):
    t_["p_holm"] = float(q_)

def cld(groups, diff, order):
    """Compact letter display, insert-and-absorb (Piepho 2004)."""
    sets = [set(groups)]
    for a_, b_ in diff:
        new = []
        for S_ in sets:
            new += [S_ - {a_}, S_ - {b_}] if (a_ in S_ and b_ in S_) else [S_]
        new = [S_ for S_ in new if S_]
        sets = [S_ for i_, S_ in enumerate(new)
                if not any(S_ < T_ or (S_ == T_ and j_ < i_) for j_, T_ in enumerate(new))]
    rank = {g: r for r, g in enumerate(order)}
    # letters in order of the groups they cover; ties broken by the next group,
    # which is the naming multcompView uses
    sets.sort(key=lambda S_: sorted(rank[g] for g in S_))
    return {g: "".join("abcdefghij"[i_] for i_, S_ in enumerate(sets) if g in S_)
            for g in groups}

g_letters = cld(GK, [(t["a"], t["b"]) for t in g_tests if t["p_holm"] < 0.05], GK)
print(f"panel G letters (label permutation x {N_PERM_G}, Holm p < 0.05):")
for k, m_ in zip(GK, obs_med):
    print(f"   {k.replace(chr(10), ' '):18s} n={len(gg[k]):3d}  median {m_:.3f}  {g_letters[k]}")
for t_ in g_tests:
    print(f"      {t_['a'].replace(chr(10),' '):18s} vs {t_['b'].replace(chr(10),' '):18s}"
          f"  p={t_['p']:.4f}  p_holm={t_['p_holm']:.4f}")

# ---------- panel H --------------------------------------------------------
pairs = list(itertools.combinations(SITES, 2))
hh = {}
for s1, s2 in pairs:
    for t_ in "PW":
        hh[f"{s1}-{s2}|{t_}"] = [BC[i, j] for i, j in itertools.combinations(range(n), 2)
                                 if typ[i] == typ[j] == t_ and {site[i], site[j]} == {s1, s2}]

json.dump(dict(source=SOURCE, sample=sm, site=site.tolist(), type=typ.tolist(), hab=hab.tolist(),
               shannon=H.tolist(), alpha=alpha,
               pcoa=dict(x=CO[:, 0].tolist(), y=CO[:, 1].tolist(),
                         exp=[float(EXP[0]), float(EXP[1])]),
               permanova=dict(site=pf_site, type=pf_type, interaction=ph_int),
               panelG={k: v for k, v in gg.items()},
               panelG_letters=g_letters, panelG_tests=g_tests, bc_input=BC_INPUT,
               panelH=hh, pairs=["-".join(p) for p in pairs]),
          open("output/fig3.json", "w"), indent=1)
print("saved")
