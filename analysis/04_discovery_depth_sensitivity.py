"""
Sensitivity of the OTU results to the OTU-discovery subsampling depth.

Usage: python analysis/04_discovery_depth_sensitivity.py TABLE_A TABLE_B [LABEL_A LABEL_B]

Each TABLE is an otu_table.tsv from preprocessing/03_otu97_pipeline.sh run with a
different BUILD_N (reads per library used to build the OTUs). Every other
parameter must be identical, so the only variable is discovery depth.

Presence rules (identical to Fig. 4):
    detected in a library   >= MIN_READS (8) mapped reads
    present at a station    detected in >= 2 of that station's 3 replicates
    present in a sample type present at >= 1 station of that type

Reported per table, for plastic (P) and water (W):
    total OTUs, site-specific OTUs, OTUs at all four stations,
    plastic-only / shared / water-only, and mean richness per library
    (observed, and rarefied to the shallowest library).
"""
import sys
import numpy as np
import pandas as pd

MIN_READS, MINPREV, N_RAR, SEED = 8, 2, 10, 1
SITES = ["R1", "R2", "S1", "Sc"]


def summarise(path, label):
    raw = pd.read_csv(path, sep="\t", index_col=0)
    assert raw.shape[1] == 24, f"{path}: expected 24 libraries, found {raw.shape[1]}"
    hit = raw >= MIN_READS

    def at_station(site, typ):
        cols = [f"{site}{typ}{i}" for i in (1, 2, 3)]
        return hit[cols].sum(axis=1) >= MINPREV

    out = {"label": label, "n_otus_table": len(raw)}
    pres = {}
    for typ in "PW":
        per_site = pd.DataFrame({s: at_station(s, typ) for s in SITES})
        anywhere = per_site.any(axis=1)
        pres[typ] = anywhere
        n_total = int(anywhere.sum())
        n_single = int((per_site.sum(axis=1) == 1).sum())
        n_all4 = int(per_site.all(axis=1).sum())
        out[f"{typ}_total"] = n_total
        out[f"{typ}_site_specific"] = n_single
        out[f"{typ}_site_specific_pct"] = 100 * n_single / n_total
        out[f"{typ}_all_four"] = n_all4

    both = pres["P"] & pres["W"]
    out["plastic_only"] = int((pres["P"] & ~pres["W"]).sum())
    out["shared"] = int(both.sum())
    out["water_only"] = int((pres["W"] & ~pres["P"]).sum())
    out["union"] = out["plastic_only"] + out["shared"] + out["water_only"]

    # richness per library, observed and rarefied to the shallowest library
    rng = np.random.default_rng(SEED)
    depth = int(raw.sum().min())
    obs = hit.sum()
    rar = np.mean([[(rng.multivariate_hypergeometric(
        raw[c].to_numpy(np.int64), depth) >= MIN_READS).sum()
        for c in raw.columns] for _ in range(N_RAR)], axis=0)
    rar = pd.Series(rar, index=raw.columns)
    typ_of = np.array([c[2] for c in raw.columns])
    out["rarefy_depth"] = depth
    for typ, nm in (("P", "plastic"), ("W", "water")):
        m = typ_of == typ
        out[f"richness_obs_{nm}"] = float(obs[m].mean())
        out[f"richness_rar_{nm}"] = float(rar[m].mean())
    return out


rows = [summarise(sys.argv[1], sys.argv[3] if len(sys.argv) > 3 else "A"),
        summarise(sys.argv[2], sys.argv[4] if len(sys.argv) > 4 else "B")]

keys = [("n_otus_table", "OTUs in table", "{:,.0f}"),
        ("P_total", "Plastic OTUs (present at >=1 station)", "{:,.0f}"),
        ("W_total", "Water OTUs (present at >=1 station)", "{:,.0f}"),
        ("P_site_specific_pct", "Plastic OTUs at one station only (%)", "{:.1f}"),
        ("W_site_specific_pct", "Water OTUs at one station only (%)", "{:.1f}"),
        ("P_all_four", "Plastic OTUs at all four stations", "{:,.0f}"),
        ("W_all_four", "Water OTUs at all four stations", "{:,.0f}"),
        ("plastic_only", "Plastic-only OTUs", "{:,.0f}"),
        ("shared", "Shared OTUs", "{:,.0f}"),
        ("water_only", "Water-only OTUs", "{:,.0f}"),
        ("union", "Total (plastic-only + shared + water-only)", "{:,.0f}"),
        ("richness_obs_plastic", "Mean observed richness, plastic", "{:,.1f}"),
        ("richness_obs_water", "Mean observed richness, water", "{:,.1f}"),
        ("richness_rar_plastic", "Mean rarefied richness, plastic", "{:,.1f}"),
        ("richness_rar_water", "Mean rarefied richness, water", "{:,.1f}")]

w = max(len(lbl) for _, lbl, _ in keys)
print(f"{'':{w}}  {rows[0]['label']:>12}  {rows[1]['label']:>12}")
for k, lbl, fmt in keys:
    print(f"{lbl:{w}}  {fmt.format(rows[0][k]):>12}  {fmt.format(rows[1][k]):>12}")
print(f"\nrarefaction depth: {rows[0]['rarefy_depth']:,} and "
      f"{rows[1]['rarefy_depth']:,} reads")
