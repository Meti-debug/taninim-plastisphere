"""EPI2ME wf-metagenomics species counts -> data/counts.csv.

    python preprocessing/02_epi2me_counts_to_table.py EPI2ME_CSV BARCODE_MAP

INPUT   EPI2ME_CSV    wf-metagenomics-counts-species.csv (Kraken2/Bracken):
                      species, barcode01..barcode24, total, rank columns, tax
        BARCODE_MAP   two tab-separated columns, barcode and library name in
                      the CURRENT scheme: barcode04<TAB>R1P1
OUTPUT  data/counts.csv   tax, genus, species, then one column per library

LIBRARY LABELS
    The table carries the station labels used during sequencing, one position
    upstream of the labels used in the manuscript: R2 for R1, R3 for R2, and
    S1P4 for S1P1. This script writes those labels unchanged, so the table it
    produces matches the deposited Supplementary Data, and the analysis
    scripts convert them on read. The mapping is given in data/README.md.
"""
import sys
import pandas as pd

src, bmap = sys.argv[1:3]
k = pd.read_csv(src)
m = pd.read_csv(bmap, sep="\t", header=None, names=["barcode", "library"])
m["barcode"] = m.barcode.str.strip(); m["library"] = m.library.str.strip()
assert m.barcode.is_unique and m.library.is_unique, "duplicate barcode or library in map"

missing = sorted(set(m.barcode) - set(k.columns))
assert not missing, f"barcodes in the map but not in the EPI2ME table: {missing}"
extra = sorted(c for c in k.columns if c.startswith("barcode") and c not in set(m.barcode))
if extra:
    print("note: barcodes in the EPI2ME table with no library, left out:", ", ".join(extra))

def legacy(lib):
    """current label -> label used in the deposited Kraken2 table"""
    st = {"R1": "R2", "R2": "R3"}.get(lib[:2], lib[:2])
    out = st + lib[2:]
    return "S1P4" if out == "S1P1" else out

ren = {b: legacy(l) for b, l in zip(m.barcode, m.library)}
out = k[["tax", "genus", "species"] + list(m.barcode)].rename(columns=ren)
libs = sorted(ren.values())
out = out[["tax", "genus", "species"] + libs]
out[libs] = out[libs].fillna(0).astype(int)
out.to_csv("data/counts.csv", index=False)

print(f"{len(out):,} species x {len(libs)} libraries -> data/counts.csv")
print("reads per library:")
for b, l in zip(m.barcode, m.library):
    print(f"  {b}  {l:5s} (written as {ren[b]:5s})  {int(out[ren[b]].sum()):>9,}")
