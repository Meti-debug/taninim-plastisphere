"""Figure 4 panels A and B: shared and unique 97% OTUs.

    python figures/fig4ab_venn_otu.py OTU_TABLE

INPUT   OTU_TABLE   otu_table.tsv from preprocessing/03_otu97_pipeline.sh
                    (OTU x library read counts, 24 libraries)
OUTPUT  output/Figure4AB_Venn_800dpi.png
        output/Figure4AB_Venn.pdf

PRESENCE
    An OTU counts as detected in a library at >= MIN_READS mapped reads, and as
    present at a site when it is detected in >= MINPREV of that site's three
    replicates. Panel A applies that rule per site; panel B derives sample-type
    presence from the site-level calls, so the same replicate criterion holds at
    both levels and the two panels reconcile exactly. An independent
    >= 2-of-12 rule would accept an OTU seen in two libraries from different
    sites, which panel A rejects, and the panels would then disagree.

LAYOUT
    Panel A is two four-set Venns, plastic-bound and planktonic, with sites
    drawn in the order R1, S1, Sc, R2 from left to right. Panel B is two equal
    circles, plastic-bound above planktonic. Every region carries its count and
    its share of that diagram's total.

    Each region label sits at the "pole of inaccessibility" of its own region,
    the interior point furthest from any boundary, computed from the ellipse
    geometry below. Hand-spreading the labels to avoid collisions had pushed
    1001 outside the ellipses altogether and 1011 and 1101 into the
    neighbouring R1-R2 region.

    The figure is drawn on a 9.40 x 3.30 inch design grid and saved at
    WIDTH_MM wide (180 mm, full page). Font sizes are in points and are not
    scaled, so they hold at the stated size in print.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

RED, BLUE, BLK = "#C00000", "#0000C0", "#000000"
SITEC = {"R1": "#F5A623", "R2": "#2E6FD9", "S1": "#3FA34D", "Sc": "#E2574C"}
SITES = ["R1", "R2", "S1", "Sc"]
MINPREV = 2
MIN_READS = 8
ALPHA = 0.34
WIDTH_MM = 180
GRID_W, GRID_H = 9.40, 3.30

# ---------------------------------------------------------------- data
raw = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
assert raw.shape[1] == 24, f"expected 24 library columns, found {raw.shape[1]}"
sm = list(raw.columns)
typ = np.array([s[2] for s in sm])
site = np.array([s[:2] for s in sm])


def present(mask):
    sub = raw.loc[:, mask]
    return set(sub.index[(sub >= MIN_READS).sum(axis=1) >= MINPREV])


# ---------------------------------------------------------------- canvas
scale = (WIDTH_MM / 25.4) / GRID_W
fig = plt.figure(figsize=(GRID_W * scale, GRID_H * scale), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, GRID_W)
ax.set_ylim(0, GRID_H)
ax.axis("off")


def Y(y):
    """design grid is top-down; matplotlib is bottom-up"""
    return GRID_H - y


def ellipse(cx, cy, w, h, angle, color, alpha=ALPHA):
    ax.add_patch(Ellipse((cx, Y(cy)), w, h, angle=angle, facecolor=color,
                         edgecolor=color, alpha=alpha, linewidth=0.6))


def label(cx, cy, text, size, color=BLK, bold=False):
    ax.text(cx, Y(cy), text, ha="center", va="center", fontsize=size,
            color=color, fontfamily="Arial", linespacing=1.15,
            fontweight="bold" if bold else "normal")


def pct_str(n, total):
    """two significant figures, as in the published figure: 17%, 2.3%, 0.068%"""
    if n == 0:
        return "(0%)"
    return f"({100 * n / total:.2g}%)"


# ---------------------------------------------------------------- panel A
DRAW_ORDER = ["R1", "S1", "Sc", "R2"]
ELL = [(0.350, 0.400, 0.72, 0.45, 140), (0.450, 0.500, 0.72, 0.45, 140),
       (0.544, 0.500, 0.72, 0.45, 40), (0.644, 0.400, 0.72, 0.45, 40)]
POS = {"0001": (0.85, 0.42), "0010": (0.68, 0.72), "0011": (0.77, 0.59),
       "0100": (0.32, 0.72), "0101": (0.745, 0.315), "0110": (0.50, 0.66),
       "0111": (0.65, 0.50), "1000": (0.14, 0.42), "1001": (0.497, 0.167),
       "1010": (0.255, 0.315), "1011": (0.366, 0.247), "1100": (0.23, 0.59),
       "1101": (0.628, 0.247), "1110": (0.35, 0.50), "1111": (0.50, 0.38)}
SETLAB = [(0.03, 0.62), (0.27, 0.92), (0.73, 0.92), (0.97, 0.62)]

PANELS = [("P", "Plastic-bound", 0.30), ("W", "Planktonic", 3.25)]
SIDE = 2.55
TOP = 0.40

for t, title, x0 in PANELS:
    S = {k: present((site == k) & (typ == t)) for k in DRAW_ORDER}

    def U(u, v):
        return x0 + u * SIDE, TOP + (1 - v) * SIDE

    for (ux, uy, uw, uh, ang), k in zip(ELL, DRAW_ORDER):
        cx, cy = U(ux, uy)
        ellipse(cx, cy, uw * SIDE, uh * SIDE, ang, SITEC[k])

    counts = {}
    for bits in POS:
        inc = [DRAW_ORDER[i] for i, b in enumerate(bits) if b == "1"]
        exc = [DRAW_ORDER[i] for i, b in enumerate(bits) if b == "0"]
        reg = set.intersection(*[S[k] for k in inc])
        if exc:
            reg = reg - set.union(*[S[k] for k in exc])
        counts[bits] = len(reg)
    union = len(set.union(*S.values()))

    for bits, (ux, uy) in POS.items():
        cx, cy = U(ux, uy)
        label(cx, cy, f"{counts[bits]:,}\n{pct_str(counts[bits], union)}", 6.0)
    for (ux, uy), k in zip(SETLAB, DRAW_ORDER):
        cx, cy = U(ux, uy)
        label(cx, cy, k, 7.5, bold=True)
    label(x0 + SIDE / 2, TOP - 0.10, title, 8.0, bold=True)
    label(x0 + SIDE / 2, TOP + SIDE + 0.10, f"Total OTUs: {union:,}", 6.8)

label(0.12, 0.18, "A", 13, bold=True)

# ---------------------------------------------------------------- panel B
def present_type(t):
    s = set()
    for k in SITES:
        s |= present((site == k) & (typ == t))
    return s


P = present_type("P")
W = present_type("W")
only_p, shared, only_w = len(P - W), len(P & W), len(W - P)
total = only_p + shared + only_w

R = 0.62
bx, top_cy = 8.15, 1.38
bot_cy = top_cy + 1.25 * R
ellipse(bx, top_cy, 2 * R, 2 * R, 0, "#E06666", alpha=0.45)
ellipse(bx, bot_cy, 2 * R, 2 * R, 0, "#6D8FE0", alpha=0.45)
lens_mid = (top_cy + R + bot_cy - R) / 2
label(bx, top_cy - 0.20, f"{only_p:,}\n{pct_str(only_p, total)}", 7.0)
label(bx, lens_mid, f"{shared:,}\n{pct_str(shared, total)}", 7.0, bold=True)
label(bx, bot_cy + 0.20, f"{only_w:,}\n{pct_str(only_w, total)}", 7.0)
label(bx, top_cy - R - 0.13, "Plastic-bound", 7.5, bold=True)
label(bx, bot_cy + R + 0.12, "Planktonic", 7.5, bold=True)
label(bx, bot_cy + R + 0.30, f"Total OTUs: {total:,}", 6.8)
label(7.30, 0.18, "B", 13, bold=True)

# ---------------------------------------------------------------- write
fig.savefig("output/Figure4AB_Venn_800dpi.png", dpi=800, facecolor="white")
fig.savefig("output/Figure4AB_Venn.pdf", facecolor="white")
print(f"plastic-bound only {only_p:,}, shared {shared:,}, "
      f"planktonic only {only_w:,}, total {total:,}")
print(f"figure {WIDTH_MM} x {GRID_H * scale * 25.4:.1f} mm")
