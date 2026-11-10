# Core microbiome heatmap
#
# Heatmap of the 15 most abundant species detected in all 12 samples of a sample type.
#
# Extracted from r/figures_as_run.R, lines 179-269, which is the session code
# deposited verbatim. Paths were made relative; nothing else was changed.
#
# Run from the repository root:
#   Rscript r/04_core_microbiome_heatmap.R

source("r/00_paths.R")

#Figure 6. Core microbiome composition across plastisphere and planktonic communities
library(dplyr)
library(pheatmap)

# ── 1. Load data ─────────────────────────────────────────────
count_path <- file.path(DATA, "counts.csv")
meta_path  <- file.path(DATA, "metadata.tsv")

meta <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sp_tab <- read.csv(count_path, header = TRUE, check.names = FALSE)

# ── 2. Build count matrix ────────────────────────────────────
sample_cols         <- meta$SampleID
sp_counts           <- sp_tab[, sample_cols]
rownames(sp_counts) <- sp_tab$species
sp_counts           <- sp_counts[, meta$SampleID]

# ── 3. Rename sites ──────────────────────────────────────────
meta$Site <- dplyr::recode(meta$Site, "R2" = "R1", "R3" = "R2")

# ── 4. Core microbiome (present in all samples) ──────────────
core_mask   <- rowSums(sp_counts > 0) == ncol(sp_counts)
core_counts <- sp_counts[core_mask, ]

# ── 5. Relative abundance ────────────────────────────────────
core_rel <- sweep(core_counts, 2, colSums(sp_counts), "/") * 100

# ── 6. Order: Plastic LEFT, Water RIGHT ──────────────────────
meta_ordered <- meta %>%
  mutate(SampleType = factor(SampleType, levels = c("Plastic", "Water"))) %>%
  arrange(SampleType, Site)

core_rel_ordered <- core_rel[, meta_ordered$SampleID]

# ── 7. Top 15 most abundant core species ─────────────────────
top_n      <- 15
mean_abund <- rowMeans(core_rel_ordered)
core_top   <- core_rel_ordered[order(mean_abund, decreasing = TRUE)[1:top_n], ]

# ── 8. Column annotations ────────────────────────────────────
ann_col <- data.frame(
  SampleType = meta_ordered$SampleType,
  Site       = meta_ordered$Site,
  row.names  = meta_ordered$SampleID
)

ann_colors <- list(
  SampleType = c(Plastic = "#8B0000", Water = "#00008B"),
  Site       = c(R1 = "#E69F00", R2 = "#56B4E9", S1 = "#009E73", Sc = "#CC79A7")
)

# ── 9. Log transform ─────────────────────────────────────────
core_top[is.na(core_top)] <- 0
core_log <- log10(core_top + 0.001)

# ── 10. Clean Red → Green Palette ────────────────────────────
core_pal <- colorRampPalette(
  c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
    "#F7F7F7",
    "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")
)(100)

# ── 11. Italic species names ─────────────────────────────────
italic_labels <- sapply(rownames(core_top), function(x) {
  as.expression(bquote(italic(.(x))))
})

# ── 12. Final Heatmap with Italics ───────────────────────────
pheatmap(
  core_log,
  annotation_col      = ann_col,
  annotation_colors   = ann_colors,
  cluster_cols        = FALSE,
  cluster_rows        = TRUE,
  show_colnames       = FALSE,
  show_rownames       = TRUE,
  fontsize_row        = 10,
  fontsize            = 12,
  color               = core_pal,
  border_color        = NA,
  main                = "Core Microbiome — Plastisphere vs Planktonic",
  angle_col           = 45,
  gaps_col            = sum(meta_ordered$SampleType == "Plastic"),
  na_col              = "grey90",
  labels_row          = italic_labels   # ← Italic species names
)
 




