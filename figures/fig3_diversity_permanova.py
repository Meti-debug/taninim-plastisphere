"""Figure 3: alpha diversity, ordination, and Bray-Curtis structure.

INPUT   output/fig3.json        written by analysis/01_diversity_permanova.py
OUTPUT  output/Figure3_800dpi.png
        output/Figure3.pdf

All statistics are computed in analysis/01; this script only draws them, so
editing a number here would put the figure out of step with the table.  Change
the analysis and rerun both.

Layout note: panel sizes are set in figure-fraction coordinates at the top of
the file.  The assertions before savefig check that no text overflows the
canvas and that no label collides with the caption - they are there because
these collisions are invisible at screen resolution and obvious at 800 dpi.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, json
from matplotlib.patches import Polygon, Rectangle
from scipy.spatial import ConvexHull
from matplotlib import rcParams
rcParams["font.family"] = "Liberation Sans"
rcParams["pdf.fonttype"] = 42; rcParams["svg.fonttype"] = "none"

RED, BLUE, BLK = "#C00000", "#0000C0", "#000000"
SITEC = {"R1": "#4C9BD1", "R2": "#2FA88A", "S1": "#F0812C", "Sc": "#8B5FBF"}
GCOL = ["#A01818", "#F2A7A7", "#1C2A80", "#A9C1EC", "#555555", "#C8C8C8"]
OUT = "output/"
SITES = ["R1", "R2", "S1", "Sc"]
J = json.load(open("output/fig3.json"))
SRC = J.get("source", "Kraken2/Bracken")
rng = np.random.default_rng(3)

def stars(p): return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"
def pstr(p): return f"{p:.2e}".replace("e-0", "e-")

def frame(ax, lw=1.25):
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(True); ax.spines[s].set_linewidth(lw)
    ax.tick_params(width=0.9, length=2.8, pad=2, top=False, right=False, labelsize=7)

def boxes(ax, data, pos, cols, width=0.55, pts=True, ptsize=6, ptalpha=1.0,
          jitter=0.10, fliers=True):
    bp = ax.boxplot(data, positions=pos, widths=width, patch_artist=True,
                    showfliers=fliers, medianprops=dict(color="0.25", lw=1.0),
                    whiskerprops=dict(color=BLK, lw=0.8),
                    capprops=dict(color=BLK, lw=0.8),
                    flierprops=dict(marker="o", ms=2.0, markerfacecolor="none",
                                    markeredgecolor="0.35", markeredgewidth=0.5))
    for b, c in zip(bp["boxes"], cols):
        b.set(facecolor=c, edgecolor=BLK, linewidth=0.8)
    if pts:
        for v, x in zip(data, pos):
            v = np.asarray(v, float)
            ax.scatter(x + rng.uniform(-jitter, jitter, len(v)), v, s=ptsize,
                       color=BLK, alpha=ptalpha, zorder=5, linewidths=0)
    return bp

FW, FH = 11.2, 5.85
fig = plt.figure(figsize=(FW, FH))
XS = [0.055, 0.275, 0.495, 0.715]
PW, R1B, R2B, PH = 0.160, 0.625, 0.162, 0.290

# ---------------- A-E : Shannon -------------------------------------------
# One y-axis for all five Shannon panels: same limits, whole-number ticks, and
# the significance bracket at the same height. Scaling each panel to its own
# data made the tick steps differ between panels (1, 0.5, ...) and made
# panels look comparable when they were not.
_sh = [v for k in SITES + ["ALL"] for t_ in "PW" for v in J["alpha"][k][t_]]
SH_LO = float(np.floor(min(_sh)))
SH_BRACKET = max(_sh) + 0.35
SH_HI = SH_BRACKET + 0.75
SH_TICKS = np.arange(SH_LO, np.floor(SH_HI) + 1, 1.0)

# One y-axis for both Bray-Curtis panels (G and H).
BC_YLIM = (0.10, 1.10)
BC_TICKS = [0.2, 0.4, 0.6, 0.8, 1.0]

def shannon_panel(box, key, title, tag, ylab=False):
    ax = fig.add_axes(box)
    a, w = J["alpha"][key]["P"], J["alpha"][key]["W"]
    boxes(ax, [a, w], [1, 2], [RED, BLUE], width=0.55,
          ptsize=7 if key != "ALL" else 6)
    rg = SH_HI - SH_LO
    ax.set_ylim(SH_LO, SH_HI)
    ax.set_yticks(SH_TICKS)
    ax.set_yticklabels([f"{t:.0f}" for t in SH_TICKS])
    ax.set_xlim(0.45, 2.55)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Plastic-bound", "Planktonic"], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=7.5)
    if ylab:
        ax.set_ylabel("Shannon diversity", fontsize=8, labelpad=3)
    frame(ax)
    p = J["alpha"][key]["t_p"]
    yb = SH_BRACKET
    ax.plot([1, 1, 2, 2], [yb - 0.03 * rg, yb, yb, yb - 0.03 * rg], color=BLK, lw=0.9)
    ax.text(1.5, yb + 0.02 * rg, stars(p), ha="center", va="bottom", fontsize=8)
    ax.text(-0.32, 1.160, tag, transform=ax.transAxes, fontsize=13, fontweight="bold",
            ha="left", va="top")
    ax.set_title(title, fontsize=9, fontweight="bold", pad=4)
    return ax

for i, (k, t) in enumerate(zip(SITES, SITES)):
    shannon_panel([XS[i], R1B, PW, PH], k, t, "ABCD"[i], ylab=(i == 0))
shannon_panel([XS[0], R2B, PW, PH], "ALL", "All sites", "E", ylab=True)

# ---------------- F : PCoA -------------------------------------------------
axF = fig.add_axes([XS[1], R2B, PW, PH])
px, py = np.array(J["pcoa"]["x"]), np.array(J["pcoa"]["y"])
site = np.array(J["site"]); typ = np.array(J["type"])
# hulls group by SAMPLE TYPE, not site: plastic-bound together, planktonic together
for t_, c in [("W", BLUE), ("P", RED)]:
    pts = np.c_[px[typ == t_], py[typ == t_]]
    v = pts[ConvexHull(pts).vertices]
    axF.add_patch(Polygon(v, closed=True, facecolor=c, alpha=0.13,
                          edgecolor="none", zorder=1))
    axF.add_patch(Polygon(v, closed=True, facecolor="none", edgecolor=c,
                          alpha=0.50, linewidth=0.8, zorder=2))
for s_ in SITES:
    for t_, mk in [("P", "o"), ("W", "^")]:
        k = (site == s_) & (typ == t_)
        axF.scatter(px[k], py[k], marker=mk, s=26, facecolor=SITEC[s_],
                    edgecolor=BLK, linewidth=0.5, zorder=4)
xr = px.max() - px.min(); yr = py.max() - py.min()
axF.set_xlim(px.min() - 0.10 * xr, px.max() + 0.08 * xr)
axF.set_ylim(py.min() - 0.10 * yr, py.max() + 0.34 * yr)
axF.set_xlabel(f"PCoA1 [{J['pcoa']['exp'][0]*100:.0f}%]", fontsize=8, labelpad=2)
axF.set_ylabel(f"PCoA2 [{J['pcoa']['exp'][1]*100:.1f}%]", fontsize=8, labelpad=3)
frame(axF)
PF = J["permanova"]
axF.text(0.025, 0.975,
         f"Site: R² = {PF['site']['R2']:.2f}{stars(PF['site']['p'])}\n"
         f"Type: R² = {PF['type']['R2']:.2f}{stars(PF['type']['p'])}",
         transform=axF.transAxes, fontsize=6.2, fontweight="bold", ha="left",
         va="top", linespacing=1.6)
axF.text(-0.32, 1.160, "F", transform=axF.transAxes, fontsize=13, fontweight="bold",
         ha="left", va="top")

# ---------------- G : BC partitioned --------------------------------------
axG = fig.add_axes([XS[2], R2B, PW, PH])
GK = list(J["panelG"].keys())
boxes(axG, [J["panelG"][k] for k in GK], list(range(1, 7)), GCOL, width=0.62,
      ptsize=1.6, ptalpha=0.40, jitter=0.22, fliers=False)
axG.set_xlim(0.4, 6.6); axG.set_ylim(*BC_YLIM)
axG.set_xticks(range(1, 7))
axG.set_xticklabels(["Same site", "Different sites"] * 3, rotation=90,
                    fontweight="bold")
# letters from analysis/01: tests on per-sample means, BH q < 0.05
for x_, k in enumerate(GK, start=1):
    if "panelG_letters" in J:
        axG.text(x_, max(J["panelG"][k]) + 0.025, J["panelG_letters"][k],
                 ha="center", va="bottom", fontsize=7.5, fontweight="bold")
for x_, g_ in zip((1.5, 3.5, 5.5), ("Plastic–\nPlastic", "Water–\nWater", "Plastic–\nWater")):
    axG.text(x_, -0.37, g_, transform=axG.get_xaxis_transform(), ha="center",
             va="top", fontsize=6.0, fontweight="bold", linespacing=0.95)
axG.set_yticks(BC_TICKS)
axG.set_yticklabels([f"{t:.1f}" for t in BC_TICKS], fontsize=7)
axG.set_ylabel("Bray–Curtis dissimilarity", fontsize=8, labelpad=3)
frame(axG)
axG.tick_params(axis="x", labelsize=6.0, length=0, pad=2)
axG.text(-0.32, 1.160, "G", transform=axG.transAxes, fontsize=13, fontweight="bold",
         ha="left", va="top")

# ---------------- H : pairwise, same type ---------------------------------
axH = fig.add_axes([XS[3], R2B, PW, PH])
pairs = J["pairs"]
dat, pos, cols = [], [], []
for i, pr in enumerate(pairs):
    for off, t_, c in [(-0.19, "P", RED), (0.19, "W", BLUE)]:
        dat.append(J["panelH"][f"{pr}|{t_}"]); pos.append(i + 1 + off); cols.append(c)
boxes(axH, dat, pos, cols, width=0.33, ptsize=5)
axH.set_xlim(0.45, 6.55); axH.set_ylim(*BC_YLIM)
axH.set_xticks(range(1, 7))
axH.set_xticklabels(pairs, fontsize=7, rotation=45, ha="right", rotation_mode="anchor")
axH.set_yticks(BC_TICKS)
axH.set_yticklabels([f"{t:.1f}" for t in BC_TICKS], fontsize=7)
frame(axH)
PI = PF["interaction"]
axH.text(-0.32, 1.160, "H", transform=axH.transAxes, fontsize=13, fontweight="bold",
         ha="left", va="top")

# ---------------- legends --------------------------------------------------
lg = fig.add_axes([0.893, 0.10, 0.100, 0.82]); lg.axis("off")
lg.set_xlim(0, 1); lg.set_ylim(0, 1)
def swatch(y, col, lab, fs=7):
    lg.add_patch(Rectangle((0.02, y - 0.012), 0.085, 0.024, facecolor=col,
                           edgecolor=BLK, lw=0.6, transform=lg.transAxes, clip_on=False))
    lg.text(0.135, y, lab, transform=lg.transAxes, fontsize=fs, va="center", ha="left")

lg.text(0.02, 0.905, "Sample type", fontsize=7.5, fontweight="bold", va="center",
        transform=lg.transAxes)
swatch(0.862, RED, "Plastic-bound"); swatch(0.822, BLUE, "Planktonic")

lg.text(0.02, 0.560, "Sampling site", fontsize=7.5, fontweight="bold", va="center",
        transform=lg.transAxes)
for i, s_ in enumerate(SITES):
    swatch(0.517 - i * 0.040, SITEC[s_], s_)

for i, (mk, lab) in enumerate([("o", "Plastic-bound"), ("^", "Planktonic")]):
    lg.scatter([0.062], [0.290 - i * 0.045], marker=mk, s=26, facecolor="0.72",
               edgecolor=BLK, linewidth=0.5, transform=lg.transAxes, clip_on=False)
    lg.text(0.135, 0.290 - i * 0.045, lab, transform=lg.transAxes, fontsize=7,
            va="center", ha="left")

# Test details (Welch t, PERMANOVA settings, S = same site, D = different site)
# go in the figure caption, not on the figure.

fig.canvas.draw()
inv = fig.transFigure.inverted()
ftop = 0.0
for t in fig.texts:
    bb = inv.transform(t.get_window_extent(fig.canvas.get_renderer()))
    assert bb[0][0] > -0.001 and bb[1][0] < 1.001, ("OVERFLOW", t.get_text())
    ftop = max(ftop, bb[1][1])
for ax_ in fig.axes:
    for t_ in ax_.texts:                      # nothing may fall off the canvas
        bb = inv.transform(t_.get_window_extent(fig.canvas.get_renderer()))
        assert bb[0][1] > -0.001, ("TEXT BELOW CANVAS", t_.get_text())
    for lbl in ax_.get_xticklabels():
        if not lbl.get_text():
            continue
        bb = inv.transform(lbl.get_window_extent(fig.canvas.get_renderer()))
        assert bb[0][1] > ftop + 0.003, ("TICK / FOOTNOTE COLLISION", lbl.get_text(),
                                         round(bb[0][1], 4), round(ftop, 4))
print("layout ok")
fig.savefig(OUT + "Figure3_800dpi.png", dpi=800, facecolor="white")
fig.savefig(OUT + "Figure3.pdf", facecolor="white")
print("done")

