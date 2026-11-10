# Eskapee heatmap
#
# Heatmap of ESKAPEE species abundance by site and sample type.
#
# Extracted from r/figures_as_run.R, lines 81-178, which is the session code
# deposited verbatim. Paths were made relative; nothing else was changed.
#
# Run from the repository root:
#   Rscript r/05_eskapee_heatmap.R

source("r/00_paths.R")

#Figure 10. Distribution of ESKAPEE pathogens along the river–estuary–sea continuum. 

library(pheatmap)
library(dplyr)

# ── 1. Load data ─────────────────────────────────────────────
file_path <- file.path(DATA, "counts.csv")
sp_tab <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# ── 2. Clean count matrix ────────────────────────────────────
tax_cols <- c("kingdom", "phylum", "class", "order", "family", 
              "genus", "tax", "superkingdom", "total")

count_mat <- sp_tab %>%
  select(-any_of(tax_cols))

species_names <- count_mat[[1]]
count_mat <- count_mat[, -1]
count_mat[] <- lapply(count_mat, as.numeric)
rownames(count_mat) <- species_names

# ── 3. Extract ESKAPEE ───────────────────────────────────────
eskape_species <- c(
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Staphylococcus aureus",
  "Acinetobacter baumannii",
  "Pseudomonas aeruginosa",
  "Enterococcus faecium",
  "Enterobacter cloacae"
)

eskape_mat <- count_mat[rownames(count_mat) %in% eskape_species, , drop = FALSE]
eskape_mat <- eskape_mat[match(eskape_species, rownames(eskape_mat)), , drop = FALSE]

# ── 4. Reorder columns: Plastic LEFT | Water RIGHT ───────────
sample_names <- colnames(eskape_mat)

plastic_samples <- sample_names[grepl("P", sample_names)]
water_samples   <- sample_names[grepl("W", sample_names)]

new_order <- c(plastic_samples, water_samples)

eskape_mat <- eskape_mat[, new_order]

# ── 5. Log10 transform ───────────────────────────────────────
eskape_log <- log10(eskape_mat + 1)

# ── 6. Your Color Palette ────────────────────────────────────
your_pal <- colorRampPalette(
  c("#B0E0E6", "#7EC8E6", "#4EB3D3", "#2B8CB3", 
    "#1D7A9C", "#0F5F7A", "#084A5E", "#003D4C")
)(100)

# ── 7. Column Annotations (with new order) ───────────────────
ann_col <- data.frame(
  Site = case_when(
    grepl("^R2", new_order) ~ "R1",
    grepl("^R3", new_order) ~ "R2",
    grepl("^S1", new_order) ~ "S1",
    grepl("^Sc", new_order) ~ "Sc"
  ),
  SampleType = case_when(
    grepl("P", new_order) ~ "Plastic-bound",
    grepl("W", new_order) ~ "Planktonic"
  ),
  row.names = new_order
)

ann_colors <- list(
  Site = c("R1" = "#E31A1C", "R2" = "#1F78B4", 
           "S1" = "#33A02C", "Sc" = "#6A3D9A"),
  SampleType = c("Plastic-bound" = "#8B0000", "Planktonic" = "#00008B")
)

# ── 8. Italic species names ──────────────────────────────────
italic_labels <- sapply(eskape_species, function(x) {
  as.expression(bquote(italic(.(x))))
})

# ── 9. Final Heatmap ─────────────────────────────────────────
pheatmap(eskape_log,
         color               = your_pal,
         cluster_rows        = FALSE,
         cluster_cols        = FALSE,
         border_color        = NA,
         main                = "log10(count + 1) of ESKAPEE species",
         angle_col           = 45,
         fontsize            = 11,
         fontsize_row        = 12,
         fontsize_col        = 9,
         cellwidth           = 16,
         cellheight          = 20,
         labels_row          = italic_labels,
         annotation_col      = ann_col,
         annotation_colors   = ann_colors,
         annotation_names_col = TRUE,
         gaps_col            = length(plastic_samples))  # ← Vertical gap between Plastic & Water
