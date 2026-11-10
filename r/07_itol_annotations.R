# Itol annotations
#
# iTOL annotation files for the reference 16S tree: branch colours by site, inner ring by sample type, outer ring by family.
#
# Extracted from r/figures_as_run.R, lines 3725-4083, which is the session code
# deposited verbatim. Paths were made relative; nothing else was changed.
#
# Run from the repository root:
#   Rscript r/07_itol_annotations.R

source("r/00_paths.R")

Phylogenetic tree construction (iTOL visualization)

Taninim Plastisphere Phylogenetic Tree - iTOL Annotation Scripts
# All scripts that successfully produced the final figure
# ============================================================
# Prerequisites:
#   install.packages(c("ape", "phangorn"))
#   library(ape)
#   library(phangorn)
# ============================================================


# ── STEP 1: Load tree ────────────────────────────────────────
library(ape)
library(phangorn)

tree <- read.tree(file.path(OUT, "taninim_tree2.treefile"))
tips <- tree$tip.label

df <- data.frame(
  tip       = tips,
  sample    = sapply(strsplit(tips, "\\|"), "[", 2),
  substrate = sapply(strsplit(tips, "\\|"), "[", 3),
  site      = sapply(strsplit(tips, "\\|"), "[", 4),
  stringsAsFactors = FALSE
)


# ── STEP 2: Substrate ring (burgundy = Plastic, dark blue = Water) ──
sub_colors <- c("Plastic" = "#800020", "Water" = "#00008B")

sub_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tSubstrate",
  "COLOR\t#000000",
  "LEGEND_TITLE\tSubstrate",
  "LEGEND_SHAPES\t1\t1",
  "LEGEND_COLORS\t#800020\t#00008B",
  "LEGEND_LABELS\tPlastic\tWater",
  "STRIP_WIDTH\t80",
  "MARGIN\t10",
  "",
  "DATA"
)
writeLines(c(sub_lines, paste(df$tip, sub_colors[df$substrate], sep="\t")),
           file.path(OUT, "itol_ring_substrate.txt"))
cat("✅ output/itol_ring_substrate.txt\n")


# ── STEP 3: Site ring (4 colors) ────────────────────────────
site_colors <- c(
  "R2" = "#1565C0",   # dark blue  = River
  "R3" = "#2E7D32",   # dark green = Estuary
  "S1" = "#E65100",   # orange     = Sea
  "Sc" = "#6A1B9A"    # purple     = Control
)

site_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tSampling Site",
  "COLOR\t#000000",
  "LEGEND_TITLE\tSampling Site",
  "LEGEND_SHAPES\t1\t1\t1\t1",
  "LEGEND_COLORS\t#1565C0\t#2E7D32\t#E65100\t#6A1B9A",
  "LEGEND_LABELS\tRiver (R2)\tEstuary (R3)\tSea (S1)\tControl (Sc)",
  "STRIP_WIDTH\t60",
  "MARGIN\t10",
  "",
  "DATA"
)
writeLines(c(site_lines, paste(df$tip, site_colors[df$site], sep="\t")),
           file.path(OUT, "itol_ring_site.txt"))
cat("✅ output/itol_ring_site.txt\n")


# ── STEP 4: Branch colors by site (thick lines) ──────────────
branch_lines <- c(
  "TREE_COLORS",
  "SEPARATOR TAB",
  "",
  "DATA"
)

for (i in seq_along(tips)) {
  branch_lines <- c(branch_lines,
    paste(tips[i], "branch", site_colors[df$site[i]], "normal", "4", sep="\t")
  )
}

writeLines(branch_lines, file.path(OUT, "itol_branches_thick.txt"))
cat("✅ output/itol_branches_thick.txt\n")


# ── STEP 5: Dominant family ring ─────────────────────────────
tax <- read.csv(file.path(DATA, "counts.csv"),
                stringsAsFactors=FALSE, check.names=FALSE)

barcode_map <- c(
  "barcode01"="R2W1", "barcode02"="R2W2", "barcode03"="R2W3",
  "barcode04"="R2P1", "barcode05"="R2P2", "barcode06"="R2P3",
  "barcode07"="R3W1", "barcode08"="R3W2", "barcode09"="R3W3",
  "barcode10"="R3P1", "barcode11"="R3P2", "barcode12"="R3P3",
  "barcode13"="S1W1", "barcode14"="S1W2", "barcode15"="S1W3",
  "barcode16"="S1P4", "barcode17"="S1P2", "barcode18"="S1P3",
  "barcode19"="ScW1", "barcode20"="ScW2", "barcode21"="ScW3",
  "barcode22"="ScP1", "barcode23"="ScP2", "barcode24"="ScP3"
)
barcodes <- names(barcode_map)

dominant_family <- sapply(barcodes, function(bc) {
  family_counts <- tapply(tax[[bc]], tax$family, sum, na.rm=TRUE)
  names(which.max(family_counts))
})
names(dominant_family) <- barcode_map

family_colors <- c(
  "Comamonadaceae"     = "#1565C0",
  "Corynebacteriaceae" = "#2E7D32",
  "Paracoccaceae"      = "#E65100",
  "Erythrobacteraceae" = "#6A1B9A",
  "Arcobacteraceae"    = "#F9A825",
  "Vibrionaceae"       = "#00838F",
  "Enterobacteriaceae" = "#AD1457",
  "Roseobacteraceae"   = "#4E342E",
  "Burkholderiaceae"   = "#558B2F"
)

ring_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tDominant Family",
  "COLOR\t#000000",
  "LEGEND_TITLE\tDominant Family (per sample)",
  paste("LEGEND_SHAPES", paste(rep(1,9), collapse="\t"), sep="\t"),
  paste("LEGEND_COLORS", paste(family_colors, collapse="\t"), sep="\t"),
  paste("LEGEND_LABELS", paste(names(family_colors), collapse="\t"), sep="\t"),
  "STRIP_WIDTH\t60",
  "MARGIN\t10",
  "",
  "DATA"
)

for (i in seq_along(tips)) {
  samp  <- df$sample[i]
  fam   <- dominant_family[samp]
  color <- family_colors[fam]
  if (is.na(color)) color <- "#CCCCCC"
  ring_lines <- c(ring_lines, paste(tips[i], color, sep="\t"))
}

writeLines(ring_lines, file.path(OUT, "itol_family_ring.txt"))
cat("✅ output/itol_family_ring.txt\n")


# ── STEP 6: Site legend only (for branch color legend) ───────
legend_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tSampling Site",
  "COLOR\t#000000",
  "LEGEND_TITLE\tSampling Site",
  "LEGEND_SHAPES\t1\t1\t1\t1",
  "LEGEND_COLORS\t#1565C0\t#2E7D32\t#E65100\t#6A1B9A",
  "LEGEND_LABELS\tRiver (R2)\tEstuary (R3)\tSea (S1)\tControl (Sc)",
  "STRIP_WIDTH\t0",
  "MARGIN\t0",
  "",
  "DATA"
)
data_lines <- paste(tips, "#00000000", sep="\t")
writeLines(c(legend_lines, data_lines), file.path(OUT, "itol_site_legend.txt"))
cat("✅ output/itol_site_legend.txt\n")


# ── SUMMARY ──────────────────────────────────────────────────
cat("\n============================================================\n")
cat("Files generated — drag onto iTOL in this order:\n")
cat("  1. output/itol_branches_thick.txt   → branch colors by site\n")
cat("  2. output/itol_ring_substrate.txt   → inner ring: Plastic/Water\n")
cat("  3. output/itol_family_ring.txt      → outer ring: dominant family\n")
cat("  4. output/itol_site_legend.txt      → legend for branch colors\n")
cat("============================================================\n")


Taninim Plastisphere Phylogenetic Tree - Run R Scripts
# Run this in PowerShell to generate all iTOL annotation files
# ============================================================

# Save the R script first, then run it with Rscript
$rScript = @"
library(ape)
library(phangorn)

tree <- read.tree(file.path(OUT, "taninim_tree2.treefile"))
tips <- tree\$tip.label

df <- data.frame(
  tip       = tips,
  sample    = sapply(strsplit(tips, "\\|"), "[", 2),
  substrate = sapply(strsplit(tips, "\\|"), "[", 3),
  site      = sapply(strsplit(tips, "\\|"), "[", 4),
  stringsAsFactors = FALSE
)

# ── STEP 2: Substrate ring ───────────────────────────────────
sub_colors <- c("Plastic" = "#800020", "Water" = "#00008B")
sub_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tSubstrate",
  "COLOR\t#000000",
  "LEGEND_TITLE\tSubstrate",
  "LEGEND_SHAPES\t1\t1",
  "LEGEND_COLORS\t#800020\t#00008B",
  "LEGEND_LABELS\tPlastic\tWater",
  "STRIP_WIDTH\t80",
  "MARGIN\t10",
  "",
  "DATA"
)
writeLines(c(sub_lines, paste(df\$tip, sub_colors[df\$substrate], sep="\t")),
           file.path(OUT, "itol_ring_substrate.txt"))
cat("OK output/itol_ring_substrate.txt\n")

# ── STEP 3: Site ring ────────────────────────────────────────
site_colors <- c("R2"="#1565C0","R3"="#2E7D32","S1"="#E65100","Sc"="#6A1B9A")
site_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tSampling Site",
  "COLOR\t#000000",
  "LEGEND_TITLE\tSampling Site",
  "LEGEND_SHAPES\t1\t1\t1\t1",
  "LEGEND_COLORS\t#1565C0\t#2E7D32\t#E65100\t#6A1B9A",
  "LEGEND_LABELS\tRiver (R2)\tEstuary (R3)\tSea (S1)\tControl (Sc)",
  "STRIP_WIDTH\t60",
  "MARGIN\t10",
  "",
  "DATA"
)
writeLines(c(site_lines, paste(df\$tip, site_colors[df\$site], sep="\t")),
           file.path(OUT, "itol_ring_site.txt"))
cat("OK output/itol_ring_site.txt\n")

# ── STEP 4: Branch colors ────────────────────────────────────
branch_lines <- c("TREE_COLORS","SEPARATOR TAB","","DATA")
for (i in seq_along(tips)) {
  branch_lines <- c(branch_lines,
    paste(tips[i], "branch", site_colors[df\$site[i]], "normal", "4", sep="\t"))
}
writeLines(branch_lines, file.path(OUT, "itol_branches_thick.txt"))
cat("OK output/itol_branches_thick.txt\n")

# ── STEP 5: Dominant family ring ─────────────────────────────
tax <- read.csv(file.path(DATA, "counts.csv"),
                stringsAsFactors=FALSE, check.names=FALSE)
barcode_map <- c(
  "barcode01"="R2W1","barcode02"="R2W2","barcode03"="R2W3",
  "barcode04"="R2P1","barcode05"="R2P2","barcode06"="R2P3",
  "barcode07"="R3W1","barcode08"="R3W2","barcode09"="R3W3",
  "barcode10"="R3P1","barcode11"="R3P2","barcode12"="R3P3",
  "barcode13"="S1W1","barcode14"="S1W2","barcode15"="S1W3",
  "barcode16"="S1P4","barcode17"="S1P2","barcode18"="S1P3",
  "barcode19"="ScW1","barcode20"="ScW2","barcode21"="ScW3",
  "barcode22"="ScP1","barcode23"="ScP2","barcode24"="ScP3"
)
dominant_family <- sapply(names(barcode_map), function(bc) {
  fam_counts <- tapply(tax[[bc]], tax\$family, sum, na.rm=TRUE)
  names(which.max(fam_counts))
})
names(dominant_family) <- barcode_map

family_colors <- c(
  "Comamonadaceae"="#1565C0","Corynebacteriaceae"="#2E7D32",
  "Paracoccaceae"="#E65100","Erythrobacteraceae"="#6A1B9A",
  "Arcobacteraceae"="#F9A825","Vibrionaceae"="#00838F",
  "Enterobacteriaceae"="#AD1457","Roseobacteraceae"="#4E342E",
  "Burkholderiaceae"="#558B2F"
)
ring_lines <- c(
  "DATASET_COLORSTRIP","SEPARATOR TAB",
  "DATASET_LABEL\tDominant Family","COLOR\t#000000",
  "LEGEND_TITLE\tDominant Family (per sample)",
  paste("LEGEND_SHAPES", paste(rep(1,9), collapse="\t"), sep="\t"),
  paste("LEGEND_COLORS", paste(family_colors, collapse="\t"), sep="\t"),
  paste("LEGEND_LABELS", paste(names(family_colors), collapse="\t"), sep="\t"),
  "STRIP_WIDTH\t60","MARGIN\t10","","DATA"
)
for (i in seq_along(tips)) {
  samp  <- df\$sample[i]
  fam   <- dominant_family[samp]
  color <- family_colors[fam]
  if (is.na(color)) color <- "#CCCCCC"
  ring_lines <- c(ring_lines, paste(tips[i], color, sep="\t"))
}
writeLines(ring_lines, file.path(OUT, "itol_family_ring.txt"))
cat("OK output/itol_family_ring.txt\n")

# ── STEP 6: Site legend ──────────────────────────────────────
legend_lines <- c(
  "DATASET_COLORSTRIP","SEPARATOR TAB",
  "DATASET_LABEL\tSampling Site","COLOR\t#000000",
  "LEGEND_TITLE\tSampling Site",
  "LEGEND_SHAPES\t1\t1\t1\t1",
  "LEGEND_COLORS\t#1565C0\t#2E7D32\t#E65100\t#6A1B9A",
  "LEGEND_LABELS\tRiver (R2)\tEstuary (R3)\tSea (S1)\tControl (Sc)",
  "STRIP_WIDTH\t0","MARGIN\t0","","DATA"
)
writeLines(c(legend_lines, paste(tips, "#00000000", sep="\t")),
           file.path(OUT, "itol_site_legend.txt"))
cat("OK output/itol_site_legend.txt\n")

cat("\nAll files generated successfully!\n")
cat("Drag onto iTOL in this order:\n")
cat("  1. itol_branches_thick.txt\n")
cat("  2. itol_ring_substrate.txt\n")
cat("  3. itol_family_ring.txt\n")
cat("  4. itol_site_legend.txt\n")
"@

# Save R script to temp file
$rScript | Out-File -FilePath file.path(OUT, "D:\itol_generate.R") -Encoding UTF8

# Run with Rscript
& file.path(OUT, "C:\Program Files\R\R-4.4.2\bin\Rscript.exe") file.path(OUT, "D:\itol_generate.R")

#(The reason we used PowerShell/command line was specifically for IQ-TREE2 — because it's a standalone executable that needs to be called from the command line, not from R.)

### setting used in iTOL
##Tree Display Settings (Basic tab)
Tree layout: Circular
Line width: 0px (to remove internal radial guide lines)
Branch lengths: Hidden
Internal tree scale: Hidden
Tree scale box: Hidden
Node IDs: Hidden
Internal node symbols: Hidden
Leaf node symbols: Hidden
Arc gap: Removed (set to 0)
Datasets Applied (in order)
itol_branches_thick.txt → TREE_COLORS format, branch lines colored by site, thickness 4px
itol_ring_substrate.txt → DATASET_COLORSTRIP, strip width 80px, margin 10px, burgundy (#800020) = Plastic, dark blue (#00008B) = Water
itol_family_ring.txt → DATASET_COLORSTRIP, strip width 60px, margin 10px, 9 family colors
itol_site_legend.txt → DATASET_COLORSTRIP, strip width 0px (invisible, legend only)
Substrate Dataset Settings (gear icon)
Type: Color strip
Strip width: 80px
Left margin: 50px (adjusted to push symbols away from branches)
Display dataset label: No
Circular dataset label: Yes
Cover Setting (Colored ranges panel)
Cover: Full (to fill clade sectors with solid color)
Export
Format: SVG or PDF for publication quality
This is everything applied in sequence to produce the final figure!
Sonnet 4.6

