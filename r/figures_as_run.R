#Figure 2. CHROMagar colony morphology, total growth, and ESBL confirmation across sampled rivers. (A)  CHROMagar ESBL Stacked bar chart showing the proportional distribution of colony colors for each sampled river. (B) CHROMagar Ampicillin Stacked bar chart showing the proportional distribution of colony colors on medium. Numbers inside segments indicate raw colony counts where n > 10. Colony colors correspond to putative bacterial identities based on CHROMagar manufacturer guidelines: metallic blue and turquoise colonies are indicative of Escherichia coli and related Enterobacteriaceae; pink colonies suggest Klebsiella, Enterobacter, or Citrobacter spp.; purple colonies are associated with other Gram-negative organisms; yellow tinge colonies indicate non-fermenting Gram-negatives. (C) Total colony counts comparing growth on CHROMagar (total bacteria), CHROMagar Ampicillin (ampicillin-resistant bacteria), and CHROMagar ESBL (ESBL-producing bacteria) across sampling locations. (D) ESBL confirmation from river water isolates. Oxidase-positive colonies (teal) were excluded as non-Gram-negative Enterobacteriaceae. Of the remaining oxidase-negative candidates, confirmed ESBL-positive (red) and ESBL-negative (blue) isolates are shown. Taninim River yielded the highest proportion of confirmed ESBL-producing Gram-negative bacteria (24 of 25 oxidase-negative candidates), followed by Hadera River (6 of 9), while Yarkon and Alexander rivers yielded no confirmed ESBL producers.

library(ggplot2); library(dplyr); library(ggpubr)

rivers_order <- c("Yarkon River", "Alexander River", "Hadera River", "Taninim River")

amp_data <- data.frame(
  River = rep(rivers_order, each = 10),
  Color = rep(c("Pink","Turquoise","Metallic Blue","Clear","Purple",
                "White","Yellow Tinge","Purple & White Halo",
                "Dark Purple","Offwhite"), 4),
  Count = c(30,50,80,80,50,30,20,50,50,0,
            20,100,100,80,100,50,30,0,100,0,
            30,100,100,80,100,50,50,5,50,0,
            0,100,200,200,100,200,0,0,200,0))

esbl_data <- data.frame(
  River = rep(rivers_order, each = 8),
  Color = rep(c("Pink","Turquoise","Metallic Blue","Clear",
                "Purple","White","Yellow Tinge","Offwhite"), 4),
  Count = c(1,0,0,33,0,5,0,0,
            0,1,30,0,1,5,0,10,
            100,5,50,100,80,10,0,0,
            100,100,200,100,0,100,200,0))

colony_colors <- c(
  "Pink"="#FF69B4","Turquoise"="#00BFFF","Metallic Blue"="#2F4F7F",
  "Clear"="#E8E8D0","Purple"="#8B008B","White"="#F0F0F0",
  "Yellow Tinge"="#FFEA00","Purple & White Halo"="#BA55D3",
  "Dark Purple"="#4B0082","Offwhite"="#FAF0E6")

shared_legend <- c(
  "Pink"                = "Pink — Escherichia coli",
  "Turquoise"           = "Turquoise — Enterococcus / Serratia spp.",
  "Metallic Blue"       = "Metallic Blue — Klebsiella / Enterobacter / Citrobacter",
  "Clear"               = "Clear — Pseudomonas / Stenotrophomonas",
  "Purple"              = "Purple — Klebsiella / Enterobacter / Citrobacter",
  "White"               = "White — Acinetobacter / Others",
  "Yellow Tinge"        = "Yellow Tinge — Pseudomonas aeruginosa",
  "Purple & White Halo" = "Purple & White Halo — Citrobacter / Enterobacter",
  "Dark Purple"         = "Dark Purple — Other Gram-negative",
  "Offwhite"            = "Offwhite — Other")

make_stacked <- function(data) {
  data$Color <- factor(data$Color, levels = names(colony_colors))
  ggplot(data %>% filter(Count > 0),
         aes(x = factor(River, levels = rivers_order), y = Count, fill = Color)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.45) +
    scale_fill_manual(values = colony_colors, labels = shared_legend,
                      name = NULL, drop = FALSE) +
    labs(x = NULL, y = "Number of Colonies") +
    theme_classic(base_size = 13) +
    theme(
      text = element_text(color = "black"),
      axis.text.x = element_text(size = 11.5, face = "bold", angle = 30,
                                 hjust = 1, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.y = element_text(face = "bold", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
      axis.line = element_line(color = "black", linewidth = 1.0),
      legend.position = "right",
      legend.text = element_text(size = 8.8, face = "bold.italic", color = "black"))
}

pESBL <- make_stacked(esbl_data)   # now panel A
pAmp  <- make_stacked(amp_data)    # now panel B

final <- ggarrange(pESBL, pAmp, ncol = 2, nrow = 1,
                   labels = c("A", "B"),
                   common.legend = TRUE, legend = "right",
                   font.label = list(size = 16, face = "bold", color = "black"))

ggsave("Figure2_ESBL_Amp_Final.png", final,
       width = 14, height = 6.5, dpi = 600, bg = "white")
ggsave("Figure2_ESBL_Amp_Final.pdf", final,
       width = 14, height = 6.5, device = cairo_pdf, bg = "white")

dim(png::readPNG("Figure2_ESBL_Amp_Final.png"))[2:1]


#Figure 10. Distribution of ESKAPEE pathogens along the river–estuary–sea continuum. 

library(pheatmap)
library(dplyr)

# ── 1. Load data ─────────────────────────────────────────────
file_path <- "C:/Users/metig/Downloads/Lig_wf-metagenomics-counts-species.csv"
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
#Figure 6. Core microbiome composition across plastisphere and planktonic communities
library(dplyr)
library(pheatmap)

# ── 1. Load data ─────────────────────────────────────────────
count_path <- "C:/Users/metig/Downloads/Lig_wf-metagenomics-counts-species.csv"
meta_path  <- "C:/Users/metig/Downloads/Ligation_all_metadata (1).txt"

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
 




#Taninim River as a Hotspot

library(ggplot2)
library(dplyr)
library(ggpubr)

# ── Colors ────────────────────────────────────────────────────────────────────
esbl_colors <- c(
  "Pink"         = "#FF69B4",
  "Metallic blue"= "#2F4F7F",
  "Clear"        = "#FAF0E6",
  "Purple"       = "#6A0DAD",
  "Turquoise"    = "#00BFFF",
  "White"        = "#FFFFFF",
  "Yellow Tinge" = "#FFFF00",
  "Offwhite"     = "#FFF5DC"
)

amp_colors <- c(
  "Pink"               = "#FF69B4",
  "Turquoise"          = "#00BFFF",
  "Metallic Blue"      = "#2F4F7F",
  "Clear"              = "#FAF0E6",
  "Purple"             = "#8B008B",
  "White"              = "#FFFFFF",
  "Yellow Tinge"       = "#FFFF66",
  "Purple & White Halo"= "#9966CC",
  "Dark Purple"        = "#2D0052",
  "Offwhite"           = "#FFF5DC"
)

total_colors <- c(
  "Total Chromagar"  = "#4472C4",
  "Total Ampicillin" = "#ED7D31",
  "Total ESBL"       = "#A5A5A5"
)

sites <- c("Yarkon", "Alexander River",
           "Hadera River", "Taninim River")

# ── Helper ────────────────────────────────────────────────────────────────────
to_pct <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (s == 0) return(rep(0, length(x)))
  round(x / s * 100, 1)
}
avg <- function(...) mean(c(...), na.rm = TRUE)

# ── DATA — ESBL raw counts ────────────────────────────────────────────────────
yarkon_esbl_raw  <- c(avg(0,0,0), avg(0,0,0), avg(80,5,15),
                      avg(0,0,0), avg(0,0,0), avg(5,0,5),
                      avg(0,0,0), avg(0,0,0))
alex_esbl_raw    <- c(0,  1, 30,  0,  1,  5,  0, 10)
hadera_esbl_raw  <- c(100, 50, 100, 80,  5, 10,  0,  1)
taninim_esbl_raw <- c(100, 200, 100,  0, 100, 100, 200,  0)

esbl_counts <- data.frame(
  Site  = rep(sites, each = 8),
  Color = rep(c("Pink","Metallic blue","Clear","Purple",
                "Turquoise","White","Yellow Tinge","Offwhite"), 4),
  Count = c(yarkon_esbl_raw, alex_esbl_raw,
            hadera_esbl_raw, taninim_esbl_raw)
)

esbl_df <- esbl_counts %>%
  group_by(Site) %>%
  mutate(Pct = round(Count / sum(Count) * 100, 1)) %>%
  ungroup()

esbl_df$Site  <- factor(esbl_df$Site,  levels = sites)
esbl_df$Color <- factor(esbl_df$Color,
  levels = c("Pink","Metallic blue","Clear","Purple",
             "Turquoise","White","Yellow Tinge","Offwhite"))

# ── DATA — AMP raw counts ─────────────────────────────────────────────────────
yarkon_amp_raw  <- c(avg(10,80,0),  avg(30,0,70),   avg(80,100,70),
                     avg(80,100,100), avg(50,80,0),  avg(30,50,100),
                     avg(20,0,0),   avg(20,50,130),  avg(50,100,0),
                     avg(0,0,0))
alex_amp_raw    <- c(20, 100, 100,  80, 100,  50,  30,   0, 100,  0)
hadera_amp_raw  <- c(30, 100, 100,  80, 100,  50,  50,   5,  50,  0)
taninim_amp_raw <- c( 0, 100, 200, 200, 100, 200,   0,   0, 200,  0)

amp_counts <- data.frame(
  Site  = rep(sites, each = 10),
  Color = rep(c("Pink","Turquoise","Metallic Blue","Clear","Purple",
                "White","Yellow Tinge","Purple & White Halo",
                "Dark Purple","Offwhite"), 4),
  Count = c(yarkon_amp_raw, alex_amp_raw,
            hadera_amp_raw, taninim_amp_raw)
)

amp_df <- amp_counts %>%
  group_by(Site) %>%
  mutate(Pct = round(Count / sum(Count) * 100, 1)) %>%
  ungroup()

amp_df$Site  <- factor(amp_df$Site,  levels = sites)
amp_df$Color <- factor(amp_df$Color,
  levels = c("Pink","Turquoise","Metallic Blue","Clear","Purple",
             "White","Yellow Tinge","Purple & White Halo",
             "Dark Purple","Offwhite"))

# ── DATA — Total growth ───────────────────────────────────────────────────────
total_df <- data.frame(
  Site  = rep(sites, each = 3),
  Type  = rep(c("Total Chromagar","Total Ampicillin","Total ESBL"), 4),
  Value = c(
    avg(81,10,20), avg(370,560,470), avg(81,10,20),
    900,           580,              47,
    800,           560,              345,
    1000,          1000,             800
  )
)
total_df$Site <- factor(total_df$Site, levels = sites)
total_df$Type <- factor(total_df$Type,
  levels = c("Total Chromagar","Total Ampicillin","Total ESBL"))

# ── SHARED THEME ──────────────────────────────────────────────────────────────
shared_theme <- theme_classic(base_size = 11) +
  theme(
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA,
                                     linewidth = 1.2),
    axis.line         = element_blank(),
    axis.text.x       = element_text(size = 9, color = "black",
                                     face = "bold", angle = 30,
                                     hjust = 1),
    axis.text.y       = element_text(size = 9, color = "black"),
    axis.title        = element_text(size = 10, color = "black",
                                     face = "bold"),
    plot.title        = element_text(size = 11, face = "bold",
                                     hjust = 0.5),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.key.size   = unit(0.4, "cm"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.margin       = margin(8, 8, 8, 8)
  )

# ── PANEL A — ESBL with counts ────────────────────────────────────────────────
# Label color: white on dark segments, black on light segments
esbl_label_color <- c(
  "Pink"         = "black",
  "Metallic blue"= "white",
  "Clear"        = "black",
  "Purple"       = "white",
  "Turquoise"    = "black",
  "White"        = "black",
  "Yellow Tinge" = "black",
  "Offwhite"     = "black"
)

pA <- ggplot(esbl_df %>% filter(Pct > 0),
             aes(x = Site, y = Pct, fill = Color)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.3) +

  # Counts inside segments — only show if Count > 10
  geom_text(
    data = esbl_df %>% filter(Count > 10),
    aes(label = round(Count),
        color = Color),
    position = position_stack(vjust = 0.5),
    size = 2.8, fontface = "bold",
    show.legend = FALSE
  ) +

  scale_fill_manual(values  = esbl_colors,       name = "Colony color") +
  scale_color_manual(values = esbl_label_color) +

  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0),
    limits = c(0, 105)
  ) +
  labs(
    title = "Chromagar ESBL Colony Colors",
    x     = "Sampled rivers",
    y     = "% of colonies by color"
  ) +
  shared_theme

# ── PANEL B — Ampicillin with counts ─────────────────────────────────────────
amp_label_color <- c(
  "Pink"               = "black",
  "Turquoise"          = "black",
  "Metallic Blue"      = "white",
  "Clear"              = "black",
  "Purple"             = "white",
  "White"              = "black",
  "Yellow Tinge"       = "black",
  "Purple & White Halo"= "black",
  "Dark Purple"        = "white",
  "Offwhite"           = "black"
)

pB <- ggplot(amp_df %>% filter(Pct > 0),
             aes(x = Site, y = Pct, fill = Color)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.3) +

  # Counts inside segments — only show if Count > 10
  geom_text(
    data = amp_df %>% filter(Count > 10),
    aes(label = round(Count),
        color = Color),
    position = position_stack(vjust = 0.5),
    size = 2.8, fontface = "bold",
    show.legend = FALSE
  ) +

  scale_fill_manual(values  = amp_colors,        name = "Colony color") +
  scale_color_manual(values = amp_label_color) +

  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0),
    limits = c(0, 105)
  ) +
  labs(
    title = "Chromagar Ampicillin Colony Colors",
    x     = "Sampled rivers",
    y     = "% of colonies by color"
  ) +
  shared_theme

# ── PANEL C — Total growth ────────────────────────────────────────────────────
pC <- ggplot(total_df,
             aes(x = Site, y = Value, fill = Type)) +
  geom_col(
    position  = position_dodge(width = 0.65),
    width     = 0.55,
    color     = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = total_colors, name = NULL) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1100)
  ) +
  labs(
    title = "Chromagar Total Growth Comparison",
    x     = "Sampled rivers",
    y     = "Number of colonies grown"
  ) +
  shared_theme +
  theme(legend.position = "right")

# ── COMBINE ───────────────────────────────────────────────────────────────────
final <- ggarrange(pA, pB, pC,
                   ncol       = 3,
                   nrow       = 1,
                   labels     = c("A", "B", "C"),
                   font.label = list(size = 13, face = "bold"),
                   widths     = c(1, 1, 1))

print(final)

# ── SAVE ──────────────────────────────────────────────────────────────────────
ggsave(
  "C:/Users/metig/Downloads/Chromagar_combined_panels.png",
  final,
  width  = 16,
  height = 6,
  dpi    = 300,
  bg     = "white"
)

cat("Saved: C:/Users/metig/Downloads/Chromagar_combined_panels.png\n")

CHromAgar ESBL, Amp and Oxidase

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(grid)

rivers_order <- c("Yarkon River", "Alexander River", "Hadera River", "Taninim River")

# ── COLOR PALETTES ────────────────────────────────────────────────────────────
esbl_pal <- c(
  "Pink"          = "#FF6AAF",
  "Turquoise"     = "#07C0EC",
  "Metallic Blue" = "#2E4D79",
  "Clear"         = "#E8E8D0",
  "Purple"        = "#73047A",
  "White"         = "#F0F0F0",
  "Yellow Tinge"  = "#FEFA01",
  "Offwhite"      = "#D4C9A8"
)

amp_pal <- c(
  "Dark Purple"         = "#2D0A4C",
  "Purple & White Halo" = "#C590C0",
  "Yellow Tinge"        = "#FEFA01",
  "White"               = "#F0F0F0",
  "Purple"              = "#73047A",
  "Clear"               = "#E8E8D0",
  "Metallic Blue"       = "#2E4D79",
  "Turquoise"           = "#07C0EC",
  "Pink"                = "#FF6AAF",
  "Magenta"             = "#FF00FF",
  "Offwhite"            = "#D4C9A8"
)

media_pal <- c(
  "Total CHROMagar"  = "#4472C4",
  "Total Ampicillin" = "#ED7D31",
  "Total ESBL"       = "#70AD47"
)

confirm_pal <- c(
  "ESBL+ confirmed"     = "#C00000",
  "ESBL-"               = "#4472C4",
  "Oxidase+ (excluded)" = "#00B0A0"
)

# ── BASE THEME ────────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.title       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(size = 11, angle = 15, hjust = 1, color = "black"),
    axis.text.y      = element_text(size = 11, color = "black"),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.4, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.8),
    axis.line        = element_blank(),
    axis.ticks       = element_line(color = "black", linewidth = 0.8),
    plot.margin      = margin(8, 8, 8, 8),
    legend.position  = "right"
  )

# ── PANEL A: CHROMagar ESBL Colony Colors ────────────────────────────────────
esbl_df <- data.frame(
  River = rep(rivers_order, each = 8),
  Color = rep(c("Pink","Turquoise","Metallic Blue","Clear",
                "Purple","White","Yellow Tinge","Offwhite"), 4),
  Count = c(
    # Yarkon: mean of 3 replicates
    mean(c(1,0,0)), mean(c(0,0,0)), mean(c(0,0,0)),
    mean(c(80,5,15)), mean(c(0,0,0)), mean(c(0,5,5)),
    mean(c(0,0,0)), mean(c(0,0,0)),
    # Alexander
    0, 1, 1, 30, 0, 5, 0, 10,
    # Hadera
    100, 5, 50, 100, 80, 10, 0, 0,
    # Taninim
    100, 100, 200, 100, 0, 100, 200, 0
  )
) %>%
  filter(Count > 0) %>%
  mutate(
    River = factor(River, levels = rivers_order),
    Color = factor(Color, levels = rev(names(esbl_pal)))
  ) %>%
  group_by(River) %>%
  mutate(Pct = Count / sum(Count) * 100) %>%
  ungroup()

pA <- ggplot(esbl_df, aes(x = River, y = Pct, fill = Color)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = ifelse(round(Count) >= 10, round(Count), "")),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = esbl_pal, name = "Colony color", drop = FALSE) +
  scale_y_continuous(labels = percent_format(scale = 1),
                     expand = c(0, 0), limits = c(0, 101)) +
  scale_x_discrete(labels = function(x) gsub(" River", "\nRiver", x)) +
  labs(title = "A   CHROMagar ESBL Colony Colors",
       x = "Sampled rivers", y = "% of colonies by color") +
  base_theme

# ── PANEL B: CHROMagar Ampicillin Colony Colors ───────────────────────────────
amp_df <- data.frame(
  River = rep(rivers_order, each = 11),
  Color = rep(c("Dark Purple","Purple & White Halo","Yellow Tinge","White",
                "Purple","Clear","Metallic Blue","Turquoise",
                "Pink","Magenta","Offwhite"), 4),
  Count = c(
    # Yarkon: mean of 3 replicates
    mean(c(50,100,0)), mean(c(20,50,130)), mean(c(20,0,0)),
    mean(c(30,50,100)), mean(c(50,80,100)), mean(c(80,100,100)),
    mean(c(80,100,70)), mean(c(30,0,70)),  mean(c(10,80,0)),
    0, 0,
    # Alexander
    100, 0, 30, 50, 100, 80, 100, 100, 20, 0, 0,
    # Hadera
    50, 0, 50, 50, 100, 80, 100, 100, 30, 0, 0,
    # Taninim
    200, 0, 0, 200, 100, 200, 200, 100, 0, 0, 0
  )
) %>%
  filter(Count > 0) %>%
  mutate(
    River = factor(River, levels = rivers_order),
    Color = factor(Color, levels = rev(names(amp_pal)))
  ) %>%
  group_by(River) %>%
  mutate(Pct = Count / sum(Count) * 100) %>%
  ungroup()

pB <- ggplot(amp_df, aes(x = River, y = Pct, fill = Color)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3, width = 0.75) +
  geom_text(aes(label = ifelse(round(Count) >= 30, round(Count), "")),
            position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = amp_pal, name = "Colony color", drop = FALSE) +
  scale_y_continuous(labels = percent_format(scale = 1),
                     expand = c(0, 0), limits = c(0, 101)) +
  scale_x_discrete(labels = function(x) gsub(" River", "\nRiver", x)) +
  labs(title = "B   CHROMagar Ampicillin Colony Colors",
       x = "Sampled rivers", y = "% of colonies by color") +
  base_theme

# ── PANEL C: CHROMagar Total Growth Comparison ────────────────────────────────
yarkon_ca   <- c(480, 580, 400)
yarkon_amp  <- c(400, 510, 480)
yarkon_esbl <- c(80, 10, 20)

totals_df <- data.frame(
  River = factor(rep(rivers_order, each = 3), levels = rivers_order),
  Media = rep(c("Total CHROMagar","Total Ampicillin","Total ESBL"), 4),
  Mean  = c(
    mean(yarkon_ca),  mean(yarkon_amp),  mean(yarkon_esbl),
    900, 580, 47,
    800, 550, 345,
    1000, 1000, 800
  ),
  SD = c(
    sd(yarkon_ca), sd(yarkon_amp), sd(yarkon_esbl),
    0, 0, 0,
    0, 0, 0,
    0, 0, 0
  )
) %>%
  mutate(Media = factor(Media,
    levels = c("Total CHROMagar","Total Ampicillin","Total ESBL")))

pC <- ggplot(totals_df, aes(x = River, y = Mean, fill = Media)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.75, color = NA) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.8, color = "black") +
  scale_fill_manual(values = media_pal, name = NULL) +
  scale_y_continuous(limits = c(0, 1150), expand = c(0, 0),
                     breaks = seq(0, 1000, 200)) +
  scale_x_discrete(labels = function(x) gsub(" River", "\nRiver", x)) +
  labs(title = "C   CHROMagar Total Growth Comparison",
       x = "Water Source", y = "Number of Colonies Grown") +
  base_theme

# ── PANEL D: ESBL Confirmation from River Water ───────────────────────────────
confirm_df <- data.frame(
  River    = factor(rep(rivers_order, each = 3), levels = rivers_order),
  Category = rep(c("ESBL+ confirmed","ESBL-","Oxidase+ (excluded)"), 4),
  Count    = c(
    1,  2,  7,   # Yarkon
    0,  3,  3,   # Alexander
    6,  6,  3,   # Hadera
    24, 1,  0    # Taninim
  )
) %>%
  mutate(Category = factor(Category,
    levels = c("Oxidase+ (excluded)","ESBL-","ESBL+ confirmed")))

pD <- ggplot(confirm_df, aes(x = River, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = NA, width = 0.75) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_stack(vjust = 0.5),
            size = 3.2, fontface = "bold", color = "white") +
  scale_fill_manual(values = confirm_pal, name = NULL) +
  scale_y_continuous(limits = c(0, 28), expand = c(0, 0)) +
  scale_x_discrete(labels = function(x) gsub(" River", "\nRiver", x)) +
  labs(title = "D   ESBL Confirmation from River Water",
       x = "Sampled rivers", y = "Number of isolates") +
  base_theme

# ── SAVE ──────────────────────────────────────────────────────────────────────
png("Figure2_final.png", width = 14, height = 10, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  2, 2,
  widths  = unit(c(1,1), "null"),
  heights = unit(c(1,1), "null")
)))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()
cat("Done.\n")


NDMS

library(vegan)
library(ggplot2)
library(dplyr)

## 1. Abundance
abund <- read.csv(
  "~/Downloads/_updated_All_samples_abundance-genus - All_samples_abundance-genus.csv.csv",
  header = TRUE, check.names = FALSE
)
rownames(abund) <- abund$genus
sample_cols <- names(abund)[2:25]
comm <- t(as.matrix(abund[, sample_cols]))

## 2. Metadata (corrected)
meta <- read.table(
  "~/Downloads/Ligation_all_metadata (1).txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
names(meta)[names(meta) == "#SampleID"] <- "SampleID"   # now column is SampleID [file:108]

## 3. NMDS
set.seed(123)
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)  # [file:96]

site.scores <- as.data.frame(scores(nmds, display = "sites"))
site.scores$SampleID <- rownames(site.scores)  # R2W1…ScP3

plot.df <- left_join(site.scores, meta, by = "SampleID")  # joins successfully now

## 4. Plot
ggplot(plot.df,
       aes(NMDS1, NMDS2,
           color = SampleType,
           shape = Environment)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(aes(group = Environment),
               colour = "black", linewidth = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("NMDS (Bray–Curtis), stress = ",
                   signif(nmds$stress, 3)),
    color = "Sample type",
    shape = "Environment"
  )



##to adjust the color i used…

ggplot(plot.df,
       aes(NMDS1, NMDS2,
           color = Environment,   # 4 environments
           shape = SampleType)) + # Plastic / Water
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(aes(group = Environment),
               colour = "black", linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "River"   = "darkblue",
      "Estuary" = "darkorange",
      "Sea"     = "deeppink4",   # dark pink
      "Control" = "darkgreen"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", size = 1.2)
  ) +
  labs(
    title = paste0("NMDS (Bray–Curtis), stress = ",
                   signif(nmds$stress, 3)),
    color = "Environment",
    shape = "Sample type"
  )





Core Microbiome

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ── 1. Paths ──────────────────────────────────────────────────────────
count_path <- "C:\\Users\\metig\\Downloads\\Lig_wf-metagenomics-counts-species.csv"
meta_path  <- "C:\\Users\\metig\\Downloads\\Ligation_all_metadata (1).txt"

# ── 2. Read data ──────────────────────────────────────────────────────────
meta <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sp_tab <- read.csv(count_path, header = TRUE, check.names = FALSE)

# ── 3. Build count matrix ─────────────────────────────────────────────────────
sample_cols         <- meta$SampleID
sp_counts           <- sp_tab[, sample_cols]
rownames(sp_counts) <- sp_tab$species
sp_counts           <- sp_counts[, meta$SampleID]

# ── 4. Rename sites R2→R1, R3→R2 in metadata ─────────────────────────────────
meta$Site <- dplyr::recode(meta$Site,
  "R2" = "R1",
  "R3" = "R2"
)

# Verify renaming worked
cat("Sites after renaming:", unique(meta$Site), "\n")

# ── 5. Filter: core microbiome ────────────────────────────────────────────────
# Core = present in ALL samples
core_mask   <- rowSums(sp_counts > 0) == ncol(sp_counts)
core_counts <- sp_counts[core_mask, ]
cat("Core species (present in all samples):", nrow(core_counts), "\n")

# If too few, relax to 90% presence threshold instead
# core_mask <- rowSums(sp_counts > 0) >= round(0.9 * ncol(sp_counts))
# core_counts <- sp_counts[core_mask, ]

# ── 6. Relative abundance (per sample) ───────────────────────────────────────
core_rel <- sweep(core_counts, 2, colSums(sp_counts), "/") * 100

# ── 7. Order columns: Plastic LEFT, Water RIGHT ───────────────────────────────
meta_ordered <- meta %>%
  mutate(SampleType = factor(SampleType, levels = c("Plastic", "Water"))) %>%
  arrange(SampleType, Site)

core_rel_ordered <- core_rel[, meta_ordered$SampleID]

# ── 8. Keep top 15 most abundant core species ────────────────────────────────
top_n      <- 15
mean_abund <- rowMeans(core_rel_ordered)
core_top   <- core_rel_ordered[order(mean_abund, decreasing = TRUE)[1:top_n], ]

# ── 9. Column annotation bar ─────────────────────────────────────────────────
ann_col <- data.frame(
  SampleType = meta_ordered$SampleType,
  Site       = meta_ordered$Site,
  row.names  = meta_ordered$SampleID
)

# Annotation colors — note R1 and R2 now reflect renamed sites
ann_colors <- list(
  SampleType = c(Plastic = "#8B0000",
                 Water   = "#00008B"),
  Site       = c(R1 = "#E69F00",
                 R2 = "#56B4E9",
                 S1 = "#009E73",
                 Sc = "#CC79A7")
)

# ── 10. Safe log transform ────────────────────────────────────────────────────
core_top[is.na(core_top)] <- 0

core_log <- log10(core_top + 0.001)

# Remove rows where all values are identical — causes hclust to fail
all_same <- apply(core_log, 1, function(x) length(unique(x)) == 1)
if (any(all_same)) {
  cat("Removing", sum(all_same), "rows with identical values across all samples\n")
  core_log <- core_log[!all_same, ]
}

core_log <- as.matrix(core_log)

# Final check
cat("Any Inf:", any(is.infinite(core_log)), "\n")
cat("Any NaN:", any(is.nan(core_log)), "\n")
cat("Any NA: ", any(is.na(core_log)), "\n")

# ── 11. Draw heatmap ──────────────────────────────────────────────────────────
pheatmap(
  core_log,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_cols      = FALSE,        # preserve Plastic-left / Water-right order
  cluster_rows      = TRUE,         # cluster species by pattern
  show_colnames     = FALSE,        # annotation bar replaces sample IDs
  show_rownames     = TRUE,
  fontsize_row      = 9,
  fontsize          = 12,
  color             = colorRampPalette(
                        c("white", "#FEE8C8", "#E34A33", "#8B0000")
                      )(100),
  border_color      = NA,
  main              = "Core Microbiome — Plastisphere vs Planktonic",
  angle_col         = 45,
  gaps_col          = sum(meta_ordered$SampleType == "Plastic"),
  na_col            = "grey90"
)sour
ESKAPEE pathogens



library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

## 1. Read species count table
sp_tab <- read.csv("Lig_wf-metagenomics-counts-species.csv",
                   header = TRUE, check.names = FALSE)  # first col = species, then samples [file:118]

# Inspect once in R to confirm:
# head(colnames(sp_tab))

## 2. Define sample columns explicitly (edit if needed)
# Here we assume columns 2:25 are samples, rest (if any) are taxonomy.
species_col      <- 1
first_sample_col <- 2
last_sample_col  <- 25

sp_counts <- sp_tab[, first_sample_col:last_sample_col]
rownames(sp_counts) <- sp_tab[[species_col]]

## 3. ESKAPEE species list (rows to keep)
eskape_species <- c(
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Staphylococcus aureus",
  "Acinetobacter baumannii",
  "Pseudomonas aeruginosa",
  "Enterococcus faecium",
  "Enterobacter cloacae"
)

eskape_mat <- sp_counts[rownames(sp_counts) %in% eskape_species, , drop = FALSE]

# Re‑order rows to desired order in the figure
eskape_mat <- eskape_mat[match(eskape_species, rownames(eskape_mat)), , drop = FALSE]

## 4. Transform to log10(count + 1)
eskape_log <- log10(eskape_mat + 1)

## 5. Draw heatmap
pheatmap(eskape_log,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = NA,
         main = "log10(count + 1) of ESKAPEE species",
         angle_col = 45)

ggplot(vol, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(size = 1.1, alpha = 0.9) +
  geom_vline(xintercept = 0, color = "black", linewidth = 1.3) +
  scale_y_continuous(limits = c(0, 25), oob = squish) +
  scale_color_manual(values = c(
    "Plastic-enriched" = "#8B0000",
    "Water-enriched"   = "#00008B",
    "Not significant"  = "grey75"
  )) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color="black", fill=NA, linewidth=1.4),
    axis.line = element_line(color="black", linewidth=1.3),
    axis.ticks = element_line(color="black", linewidth=1.2),
    legend.title = element_blank()
  ) +
  labs(
    title = "Differential abundance (DESeq2): Plastic vs Water",
    x = "log2 fold change (Plastic / Water, shrunken)",
    y = "-log10(FDR)"
  )

ESKAPEE abundance and read counts

# ============================================================
# ESKAPEE Pathogens — Stacked Bar Chart
# Relative abundance (%) by site x substrate
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

## ---- 1) Load data ----
counts <- read_csv("C:/Users/metig/Downloads/Lig_wf-metagenomics-counts-species.csv",
                   show_col_types = FALSE)

## ---- 2) Define ESKAPEE patterns ----
eskapee_patterns <- c(
  "Enterococcus faecium",
  "Staphylococcus aureus",
  "Klebsiella pneumoniae",
  "Acinetobacter baumannii",
  "Pseudomonas aeruginosa",
  "Enterobacter",
  "Escherichia coli"
)

## ---- 3) Clean count table ----
tax_cols <- intersect(
  c("kingdom","phylum","class","order","family","genus","tax","superkingdom","total"),
  colnames(counts)
)

species_col <- colnames(counts)[1]

counts_clean <- counts %>%
  select(-all_of(tax_cols))

# Fix S1P1 -> S1P4
if ("S1P1" %in% colnames(counts_clean)) {
  counts_clean <- counts_clean %>% rename(S1P4 = S1P1)
  cat("Fixed: S1P1 -> S1P4\n")
}

## ---- 4) Sample columns ----
sample_cols <- c("R2P1","R2P2","R2P3",
                 "R2W1","R2W2","R2W3",
                 "R3P1","R3P2","R3P3",
                 "R3W1","R3W2","R3W3",
                 "S1P4","S1P2","S1P3",
                 "S1W1","S1W2","S1W3",
                 "ScP1","ScP2","ScP3",
                 "ScW1","ScW2","ScW3")
sample_cols <- sample_cols[sample_cols %in% colnames(counts_clean)]

## ---- 5) Compute per-sample total reads (from full unfiltered table) ----
sample_totals <- counts_clean %>%
  select(all_of(sample_cols)) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "SampleTotal")

## ---- 6) Filter ESKAPEE rows ----
eskapee_rows <- counts_clean %>%
  filter(str_detect(.data[[species_col]],
                    paste(eskapee_patterns, collapse = "|"))) %>%
  filter(!str_detect(.data[[species_col]], "endosymbiont"))

cat("ESKAPEE rows found:", nrow(eskapee_rows), "\n")
print(eskapee_rows[[species_col]])

## ---- 7) Assign group labels ----
eskapee_rows <- eskapee_rows %>%
  mutate(ESKAPEE_Group = case_when(
    str_detect(.data[[species_col]], "Enterococcus faecium")    ~ "E. faecium",
    str_detect(.data[[species_col]], "Staphylococcus aureus")   ~ "S. aureus",
    str_detect(.data[[species_col]], "Klebsiella pneumoniae")   ~ "K. pneumoniae",
    str_detect(.data[[species_col]], "Acinetobacter baumannii") ~ "A. baumannii",
    str_detect(.data[[species_col]], "Pseudomonas aeruginosa")  ~ "P. aeruginosa",
    str_detect(.data[[species_col]], "Escherichia coli")        ~ "E. coli",
    str_detect(.data[[species_col]], "Enterobacter")            ~ "Enterobacter spp.",
    TRUE ~ "Other"
  ))

## ---- 8) Pivot long, compute relative abundance per sample ----
long_df <- eskapee_rows %>%
  select(ESKAPEE_Group, all_of(sample_cols)) %>%
  group_by(ESKAPEE_Group) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  pivot_longer(-ESKAPEE_Group, names_to = "Sample", values_to = "Reads") %>%
  left_join(sample_totals, by = "Sample") %>%
  mutate(RelAbund = ifelse(SampleTotal > 0, (Reads / SampleTotal) * 100, 0))

## ---- 9) Add site/substrate labels ----
long_df <- long_df %>%
  mutate(
    Site = case_when(
      str_detect(Sample, "^R2") ~ "River (R2)",
      str_detect(Sample, "^R3") ~ "Estuary (R3)",
      str_detect(Sample, "^S1") ~ "Sea (S1)",
      str_detect(Sample, "^Sc") ~ "Control (Sc)"
    ),
    Substrate = case_when(
      str_detect(Sample, "P") ~ "Plastic",
      str_detect(Sample, "W") ~ "Water"
    ),
    SiteSubstrate = paste0(Site, "\n", Substrate)
  )

site_sub_order <- c(
  "River (R2)\nPlastic",   "River (R2)\nWater",
  "Estuary (R3)\nPlastic", "Estuary (R3)\nWater",
  "Sea (S1)\nPlastic",     "Sea (S1)\nWater",
  "Control (Sc)\nPlastic", "Control (Sc)\nWater"
)
long_df$SiteSubstrate <- factor(long_df$SiteSubstrate, levels = site_sub_order)

## ---- 10) Average relative abundance across replicates ----
plot_df <- long_df %>%
  group_by(ESKAPEE_Group, SiteSubstrate) %>%
  summarise(MeanRelAbund = mean(RelAbund, na.rm = TRUE), .groups = "drop")

## ---- 11) Color palette matching your legend ----
eskapee_colors <- c(
  "E. faecium"        = "#7B2D8B",
  "S. aureus"         = "#4A90D9",
  "K. pneumoniae"     = "#E91E8C",
  "A. baumannii"      = "#8B0000",
  "P. aeruginosa"     = "#6DB33F",
  "Enterobacter spp." = "#00BCD4",
  "E. coli"           = "#B0C4DE"
)

eskapee_order <- c("E. faecium", "S. aureus", "K. pneumoniae",
                   "A. baumannii", "P. aeruginosa", "Enterobacter spp.", "E. coli")

plot_df$ESKAPEE_Group <- factor(plot_df$ESKAPEE_Group, levels = eskapee_order)

## ---- 12) Plot ----
p <- ggplot(plot_df, aes(x = SiteSubstrate, y = MeanRelAbund, fill = ESKAPEE_Group)) +
  geom_bar(stat = "identity", position = "stack",
           color = "black", linewidth = 0.4, width = 0.7) +
  scale_fill_manual(
    values = eskapee_colors,
    name   = "ESKAPEE Species",
    labels = c(expression(italic("E. faecium")),
               expression(italic("S. aureus")),
               expression(italic("K. pneumoniae")),
               expression(italic("A. baumannii")),
               expression(italic("P. aeruginosa")),
               expression(italic("Enterobacter")~"spp."),
               expression(italic("E. coli")))
  ) +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  geom_vline(xintercept = c(2.5, 4.5, 6.5),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  labs(
    title    = "ESKAPEE Pathogen Relative Abundance Across Sites and Substrates",
    subtitle = "Taninim River–Mediterranean System | Mean relative abundance (%) across replicates",
    x        = NULL,
    y        = "Mean Relative Abundance (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, color = "grey40"),
    axis.text.x        = element_text(size = 10, color = "black",
                                      face = "bold", hjust = 0.5),
    axis.text.y        = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y       = element_text(size = 12, face = "bold"),
    axis.line          = element_line(color = "black", linewidth = 1.5),
    axis.ticks         = element_line(color = "black", linewidth = 1.2),
    axis.ticks.length  = unit(0.25, "cm"),
    legend.title       = element_text(face = "bold", size = 11),
    legend.text        = element_text(size = 10),
    legend.position    = "right",
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    legend.background  = element_rect(fill = "white", color = NA)
  )

print(p)
ggsave("D:/ESKAPEE_relabund_by_site_substrate.png", p,
       width = 12, height = 6, dpi = 300, bg = "white")
cat("Saved: D:/ESKAPEE_relabund_by_site_substrate.png\n")

## ---- 13) Save summary table ----
summary_table <- plot_df %>%
  pivot_wider(names_from  = SiteSubstrate,
              values_from = MeanRelAbund,
              values_fill = 0)
write.csv(summary_table, "D:/ESKAPEE_relabund_summary.csv", row.names = FALSE)
cat("Saved: D:/ESKAPEE_relabund_summary.csv\n")


### abundance

Alpha diversity per location

# Ensure all four sites are recognized
alpha_df$Site <- factor(alpha_df$Site,
                        levels = c("R2","R3","S1","Sc"))

par(mfrow = c(2, 2),     # 4 panels
    lwd   = 2,
    cex.axis = 1.1,
    cex.lab  = 1.2)

for(s in levels(alpha_df$Site)) {
  tmp <- subset(alpha_df, Site == s)
  if(nrow(tmp) < 2) next
  boxplot(Shannon ~ SampleType, data = tmp,
          main  = s,
          ylab  = "Shannon diversity",
          xlab  = "",
          col   = c("steelblue", "red"),
          border = "black",
          boxwex = 0.6,
          lwd    = 2)
}
par(mfrow = c(1, 1))

Corrected Alpha diversity
# ── Libraries ─────────────────────────────────────────────────────────────────
library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── 1. Load data ──────────────────────────────────────────────────────────────
meta <- read.table(
  "C:\\Users\\metig\\Downloads\\Ligation_all_metadata (1).txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sp_tab <- read.csv(
  "C:\\Users\\metig\\Downloads\\Lig_wf-metagenomics-counts-species.csv",
  header = TRUE, check.names = FALSE)

# ── 2. Build count matrix ─────────────────────────────────────────────────────
sp_counts <- sp_tab[, meta$SampleID]
rownames(sp_counts) <- sp_tab$species
sp_counts <- sp_counts[, meta$SampleID]
sp_counts_clean <- sp_counts[rowSums(sp_counts) > 0, ]

# ── 3. Rename sites to match manuscript ──────────────────────────────────────
meta_clean <- meta
meta_clean$Site <- dplyr::recode(meta_clean$Site,
                                  "R2" = "R1",
                                  "R3" = "R2")

cat("Sites:", table(meta_clean$Site), "\n")
cat("Types:", table(meta_clean$SampleType), "\n")

# ── 4. Calculate Shannon diversity ────────────────────────────────────────────
shannon <- diversity(t(sp_counts_clean), index = "shannon")

alpha_df <- data.frame(
  SampleID   = meta_clean$SampleID,
  SampleType = meta_clean$SampleType,
  Site       = meta_clean$Site,
  Shannon    = shannon
)

alpha_df$SampleType <- dplyr::recode(alpha_df$SampleType,
                                      "Plastic" = "Plastic-bound",
                                      "Water"   = "Planktonic")

alpha_df$SampleType <- factor(alpha_df$SampleType,
                               levels = c("Plastic-bound", "Planktonic"))

cat("Alpha df head:\n")
print(head(alpha_df))

# ── 5. Colors and theme ───────────────────────────────────────────────────────
type_cols <- c("Plastic-bound" = "#8B0000", "Planktonic" = "#00008B")

theme_pub <- theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(color = "black", size = 11,
                                 face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(color = "black", size = 11, face = "bold"),
    axis.title   = element_text(color = "black", size = 12, face = "bold"),
    axis.line    = element_line(color = "black", linewidth = 1.2),
    axis.ticks   = element_line(color = "black", linewidth = 1.0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.text  = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 11, face = "bold")
  )

# ── 6. Per-site plot function with t-test ─────────────────────────────────────
make_site_plot <- function(site_name) {
  df    <- alpha_df %>% filter(Site == site_name)
  tt    <- t.test(Shannon ~ SampleType, data = df)
  stars <- dplyr::case_when(
    tt$p.value < 0.001 ~ "***",
    tt$p.value < 0.01  ~ "**",
    tt$p.value < 0.05  ~ "*",
    TRUE               ~ "ns"
  )
  ymax  <- max(df$Shannon)
  y_seg <- ymax * 1.06
  y_txt <- ymax * 1.11

  ggplot(df, aes(x = SampleType, y = Shannon, fill = SampleType)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5,
                 width = 0.5, linewidth = 0.7) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.7,
                shape = 21, color = "black") +
    annotate("segment", x = 1, xend = 2,
             y = y_seg, yend = y_seg, linewidth = 0.6) +
    annotate("segment", x = 1, xend = 1,
             y = y_seg, yend = y_seg * 0.99, linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2,
             y = y_seg, yend = y_seg * 0.99, linewidth = 0.6) +
    annotate("text", x = 1.5, y = y_txt,
             label = paste0(stars, " (p=",
                            formatC(tt$p.value, format = "e", digits = 2), ")"),
             size = 3.5, fontface = "bold", hjust = 0.5) +
    scale_fill_manual(values = type_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(title = site_name, x = NULL, y = "Shannon diversity") +
    theme_pub +
    theme(legend.position = "none",
          axis.text.x     = element_text(angle = 45, hjust = 1),
          plot.title      = element_text(hjust = 0.5, face = "bold"))
}

p_R1 <- make_site_plot("R1")
p_R2 <- make_site_plot("R2")
p_S1 <- make_site_plot("S1")
p_Sc <- make_site_plot("Sc")

# ── 7. Overall panel with Wilcoxon ────────────────────────────────────────────
wt_all <- wilcox.test(Shannon ~ SampleType, data = alpha_df)
stars_all <- dplyr::case_when(
  wt_all$p.value < 0.001 ~ "***",
  wt_all$p.value < 0.01  ~ "**",
  wt_all$p.value < 0.05  ~ "*",
  TRUE                   ~ "ns"
)
ymax_all  <- max(alpha_df$Shannon)
y_seg_all <- ymax_all * 1.06
y_txt_all <- ymax_all * 1.11

p_overall <- ggplot(alpha_df,
                    aes(x = SampleType, y = Shannon, fill = SampleType)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5,
               width = 0.5, linewidth = 0.7) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7,
              shape = 21, color = "black") +
  annotate("segment", x = 1, xend = 2,
           y = y_seg_all, yend = y_seg_all, linewidth = 0.6) +
  annotate("segment", x = 1, xend = 1,
           y = y_seg_all, yend = y_seg_all * 0.99, linewidth = 0.6) +
  annotate("segment", x = 2, xend = 2,
           y = y_seg_all, yend = y_seg_all * 0.99, linewidth = 0.6) +
  annotate("text", x = 1.5, y = y_txt_all,
           label = paste0(stars_all, " (p=",
                          formatC(wt_all$p.value, format = "e", digits = 2), ")"),
           size = 4, fontface = "bold", hjust = 0.5) +
  scale_fill_manual(values = type_cols, name = "Sample type") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(title = "Overall Shannon diversity\nPlastic-bound vs Planktonic",
       x = NULL, y = "Shannon diversity") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5, face = "bold"))

# ── 8. Merge ──────────────────────────────────────────────────────────────────
left_panel <- (p_R1 | p_R2) / (p_S1 | p_Sc)

fig_shannon <- (left_panel | p_overall) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  ) +
  plot_layout(widths = c(2, 1))

print(fig_shannon)

ggsave("fig_shannon_diversity.pdf",
       fig_shannon, width = 14, height = 8,
       device = cairo_pdf)

ggsave("fig_shannon_diversity.png",
       fig_shannon, width = 14, height = 8, dpi = 300)

cat("Done.\n")
#Corrected Beta diversity
# ── Libraries ─────────────────────────────────────────────────────────────────
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── 1. Load data ──────────────────────────────────────────────────────────────
meta <- read.table(
  "C:\\Users\\metig\\Downloads\\Ligation_all_metadata (1).txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sp_tab <- read.csv(
  "C:\\Users\\metig\\Downloads\\Lig_wf-metagenomics-counts-species.csv",
  header = TRUE, check.names = FALSE)

# ── 2. Build count matrix ─────────────────────────────────────────────────────
sp_counts <- sp_tab[, meta$SampleID]
rownames(sp_counts) <- sp_tab$species
sp_counts <- sp_counts[, meta$SampleID]
sp_counts_clean <- sp_counts[rowSums(sp_counts) > 0, ]

# ── 3. Rename sites ───────────────────────────────────────────────────────────
meta_clean <- meta
meta_clean$Site <- dplyr::recode(meta_clean$Site,
                                  "R2" = "R1",
                                  "R3" = "R2")
cat("Sites:", table(meta_clean$Site), "\n")
cat("Types:", table(meta_clean$SampleType), "\n")

# ── 4. Bray-Curtis ────────────────────────────────────────────────────────────
bray <- vegdist(t(sp_counts_clean), method = "bray")

# ── 5. PCoA ───────────────────────────────────────────────────────────────────
pcoa_res <- cmdscale(bray, eig = TRUE, k = 2)
pcoa_var <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)

pcoa_df <- data.frame(
  PCoA1      = pcoa_res$points[, 1],
  PCoA2      = pcoa_res$points[, 2],
  SampleType = meta_clean$SampleType,
  Site       = meta_clean$Site,
  SampleID   = meta_clean$SampleID
)

# ── 6. Convex hulls ───────────────────────────────────────────────────────────
get_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]

hulls <- pcoa_df %>%
  group_by(Site) %>%
  do(get_hull(.)) %>%
  ungroup()

hull_areas <- hulls %>%
  group_by(Site) %>%
  summarise(area = abs(sum(diff(c(PCoA1, PCoA1[1])) *
                             (PCoA2 + c(PCoA2[-1], PCoA2[1])))/2)) %>%
  arrange(desc(area))

hulls <- hulls %>%
  mutate(Site = factor(Site, levels = hull_areas$Site))

# ── 7. PERMANOVA ──────────────────────────────────────────────────────────────
adonis_ST <- adonis2(bray ~ SampleType, data = meta_clean)
adonis_S  <- adonis2(bray ~ Site,       data = meta_clean)

# ── 8. Bray-Curtis long format ────────────────────────────────────────────────
bray_df <- as.data.frame(as.matrix(bray))
bray_df$Sample1 <- meta_clean$SampleID
colnames(bray_df)[1:24] <- meta_clean$SampleID

bray_long <- bray_df %>%
  pivot_longer(-Sample1,
               names_to  = "Sample2",
               values_to = "distance") %>%
  filter(Sample1 < Sample2) %>%
  mutate(
    Type1 = meta_clean$SampleType[match(Sample1, meta_clean$SampleID)],
    Site1 = meta_clean$Site[match(Sample1, meta_clean$SampleID)],
    Type2 = meta_clean$SampleType[match(Sample2, meta_clean$SampleID)],
    Site2 = meta_clean$Site[match(Sample2, meta_clean$SampleID)]
  ) %>%
  mutate(
    comparison = case_when(
      Type1 == Type2 & Site1 == Site2 ~ "Within group\n(same type, same site)",
      Type1 == Type2 & Site1 != Site2 ~ "Same type\ndiff site",
      Type1 != Type2 & Site1 == Site2 ~ "Diff type\nsame site",
      TRUE                            ~ "Diff type\ndiff site"
    ),
    site_pair = sapply(1:n(), function(i) {
      paste(sort(c(Site1[i], Site2[i])), collapse = "-")
    })
  )

# Order
pair_levels <- c("R1-R2", "R1-S1", "R1-Sc", "R2-S1", "R2-Sc", "S1-Sc")
bray_long$site_pair  <- factor(bray_long$site_pair,  levels = pair_levels)
bray_long$comparison <- factor(bray_long$comparison,
  levels = c("Within group\n(same type, same site)",
             "Same type\ndiff site",
             "Diff type\nsame site",
             "Diff type\ndiff site"))

kw_b <- kruskal.test(distance ~ comparison, data = bray_long)

cat("bray_long rows:", nrow(bray_long), "\n")

# ── 9. Colors and theme ───────────────────────────────────────────────────────
site_cols <- c(R1 = "#F28E2B", R2 = "#4E79A7", S1 = "#59A14F", Sc = "#E15759")

theme_pub <- theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(color = "black", size = 11,
                                 face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(color = "black", size = 11, face = "bold"),
    axis.title   = element_text(color = "black", size = 12, face = "bold"),
    axis.line    = element_line(color = "black", linewidth = 1.2),
    axis.ticks   = element_line(color = "black", linewidth = 1.0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.text  = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 11, face = "bold")
  )

border_theme <- theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
  axis.line    = element_blank()
)

# ── 10. Figure A — PCoA ───────────────────────────────────────────────────────
fig_beta_NMDS <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) +
  geom_polygon(data = hulls,
               aes(fill = Site, group = Site),
               color = NA, alpha = 0.25) +
  geom_point(aes(fill = Site, shape = SampleType),
             size = 4.5, stroke = 1.0, color = "black") +
  scale_shape_manual(
    values = c(Plastic = 21, Water = 24),
    labels = c(Plastic = "Plastic-bound", Water = "Planktonic")) +
  scale_fill_manual(values = site_cols,
                    breaks = c("R1","R2","S1","Sc")) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("Site: R²=", round(adonis_S$R2[1], 3),
                          ", p=", adonis_S$`Pr(>F)`[1],
                          "\nSample type: R²=", round(adonis_ST$R2[1], 3),
                          ", p=", adonis_ST$`Pr(>F)`[1]),
           hjust = -0.05, vjust = 1.4, size = 3.5, fontface = "bold") +
  labs(x     = paste0("PCoA1 [", pcoa_var[1], "%]"),
       y     = paste0("PCoA2 [", pcoa_var[2], "%]"),
       fill  = "Site", shape = "Sample type") +
  guides(
    fill  = guide_legend(override.aes = list(
      shape = 21, size  = 5,
      fill  = c("#F28E2B","#4E79A7","#59A14F","#E15759"),
      color = "black")),
    shape = guide_legend(override.aes = list(
      fill = "grey60", color = "black", size = 4))
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  border_theme

# ── 11. Figure B — Boxplot ────────────────────────────────────────────────────
fig_B <- ggplot(bray_long,
                aes(x = comparison, y = distance, fill = comparison)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5,
               width = 0.55, linewidth = 0.7) +
  scale_fill_manual(values = c(
    "Within group\n(same type, same site)" = "#2C7BB6",
    "Same type\ndiff site"                 = "#ABD9E9",
    "Diff type\nsame site"                 = "#FDAE61",
    "Diff type\ndiff site"                 = "#D7191C"
  )) +
  scale_x_discrete(labels = c(
    "Within group\n(same type, same site)" = "Within\ncommunity",
    "Same type\ndiff site"                 = "Same habitat\ndiff site",
    "Diff type\nsame site"                 = "Plastic vs Water\nsame site",
    "Diff type\ndiff site"                 = "Plastic vs Water\ndiff site"
  )) +
  annotate("text", x = 0.6, y = max(bray_long$distance) * 1.02,
           label = paste0("Kruskal-Wallis p=",
                          formatC(kw_b$p.value, format = "e", digits = 2)),
           hjust = 0, size = 3.5, fontface = "bold") +
  labs(x = NULL, y = "Bray-Curtis dissimilarity") +
  theme_pub +
  theme(legend.position = "none",
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5)) +
  border_theme

# ── 12. Figure C — Line graph ─────────────────────────────────────────────────
fig_C <- ggplot(bray_long %>% filter(!is.na(site_pair)),
                aes(x = site_pair, y = distance,
                    color = Type1, fill = Type1,
                    shape = Type1)) +
  geom_jitter(width = 0.2, size = 2.5,
              alpha = 0.5, stroke = 0.6) +
  stat_summary(fun = mean, geom = "line",
               aes(group = Type1),
               linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point",
               size = 3, shape = 21,
               color = "black") +
  scale_color_manual(values = c(Plastic = "#8B0000", Water = "#00008B"),
                     labels = c(Plastic = "Plastic-bound",
                                Water   = "Planktonic")) +
  scale_fill_manual(values  = c(Plastic = "#8B0000", Water = "#00008B"),
                    labels  = c(Plastic = "Plastic-bound",
                                Water   = "Planktonic")) +
  scale_shape_manual(values = c(Plastic = 16, Water = 17),
                     labels = c(Plastic = "Plastic-bound",
                                Water   = "Planktonic")) +
  labs(x     = "Site comparison",
       y     = "Bray-Curtis dissimilarity",
       color = "Sample type",
       fill  = "Sample type",
       shape = "Sample type") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  border_theme

# ── 13. Merge ─────────────────────────────────────────────────────────────────
fig_merged <- (fig_beta_NMDS | fig_B | fig_C) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  ) +
  plot_layout(widths = c(1.0, 0.8, 1.5))

print(fig_merged)

ggsave("fig_beta_diversity_final.pdf",
       fig_merged, width = 20, height = 6,
       device = cairo_pdf)

ggsave("fig_beta_diversity_final.png",
       fig_merged, width = 20, height = 6, dpi = 300)

cat("Done.\n")


#Alpha diversity per sample type (Plastic versus water)

library(vegan)

# Compute Shannon per sample
alpha_shannon <- diversity(t(sp_counts_clean), index = "shannon")

alpha_df <- data.frame(
  Sample     = colnames(sp_counts_clean),
  Shannon    = alpha_shannon,
  SampleType = meta_clean$SampleType,
  Site       = meta_clean$Site,
  Environment = meta_clean$Environment
)

# Plot: Water (blue) vs Plastic (red) with dark/thick lines
par(lwd = 2, cex.axis = 1.3, cex.lab = 1.4,
    col.axis = "black", col.lab = "black")

boxplot(Shannon ~ SampleType, data = alpha_df,
        ylab = "Shannon diversity",
        xlab = "",
        col    = c("steelblue", "red"),
        border = "black",
        boxwex = 0.6,
        lwd    = 2)

# Wilcoxon test
wilcox.test(Shannon ~ SampleType, data = alpha_df)


Load data
library(vegan)
library(ggplot2)
library(ggpubr)

# 1.1 Load metadata
meta_full <- read.table("Ligation_all_metadata (1).txt",
                        header = TRUE, stringsAsFactors = FALSE)

# 1.2 sp_counts: species (rows) × samples (columns), already in R
# Check sample names
colnames(sp_counts)

# 1.3 Keep only samples present in sp_counts
meta_sub <- meta_full[meta_full$SampleID %in% colnames(sp_counts), ]

# 1.4 Reorder metadata to match sp_counts column order
meta_sub <- meta_sub[match(colnames(sp_counts), meta_sub$SampleID), ]

# 1.5 Remove any rows with missing IDs or types
meta_clean <- meta_sub[!is.na(meta_sub$SampleID) & !is.na(meta_sub$SampleType), ]

# 1.6 Align count matrix to cleaned metadata
sp_counts_clean <- sp_counts[ , meta_clean$SampleID]

# Final sanity check
stopifnot(all.equal(colnames(sp_counts_clean), meta_clean$SampleID) == TRUE)

Plastic Vs water– plot + Wilcoxon


# Wilcoxon test
wilcox.test(Shannon ~ SampleType, data = alpha_df)

# Boxplot with asterisks
fig_alpha_overall <- ggboxplot(alpha_df,
                               x = "SampleType", y = "Shannon",
                               fill = "SampleType",
                               palette = c("Water" = "steelblue",
                                           "Plastic" = "red"),
                               add = "jitter") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     tip.length = 0.02,
                     size = 5) +
  ylab("Shannon diversity") +
  xlab("") +
  theme_classic(base_size = 14) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line  = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "right",
    legend.title = element_text(color = "black"),
    legend.text  = element_text(color = "black")
  )

fig_alpha_overall


Alpha diversity per site with ns

fig_alpha_site <- ggboxplot(alpha_df,
                            x = "SampleType", y = "Shannon",
                            fill = "SampleType",
                            palette = c("Water" = "steelblue",
                                        "Plastic" = "red"),
                            add = "jitter") +
  facet_wrap(~ Site, nrow = 2) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     size = 4) +
  ylab("Shannon diversity") +
  xlab("") +
  theme_classic(base_size = 13) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line  = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "right",
    legend.title = element_text(color = "black"),
    legend.text  = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black")
  )

fig_alpha_site
Beta diversity (Bray–Curtis, PERMANOVA, NMDS)
bray <- vegdist(t(sp_counts_clean), method = "bray")

PERMANOVA
# Overall SampleType effect
adonis_ST <- adonis2(bray ~ SampleType, data = meta_clean)

# Overall Site effect
adonis_S  <- adonis2(bray ~ Site, data = meta_clean)

# Full model with interaction
adonis_int <- adonis2(bray ~ SampleType * Site, data = meta_clean)
adonis_int

NMDS ordination and base‑R plot
scr <- scores(ord, display = "sites")
ord_df <- data.frame(scr,
                     SampleType = alpha_df$SampleType,
                     Site       = alpha_df$Site,
                     Environment = meta_clean$Environment,
                     stringsAsFactors = FALSE)

fig_beta_NMDS <- ggplot(ord_df, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = Environment), type = "t",
               level = 0.95, size = 0.6) +
  geom_point(aes(color = Environment,
                 shape = SampleType),
             size = 3) +
  scale_shape_manual(values = c(Water = 17, Plastic = 16)) +
  ggtitle(paste0("NMDS (Bray–Curtis), stress = ",
                 signif(ord$stress, 3))) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_classic(base_size = 13) +
  theme(
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line  = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1),
    legend.title = element_text(color = "black"),
    legend.text  = element_text(color = "black")
  )

fig_beta_NMDS

DESeq2 differential abundance, spatial gradient, SourceTracker2, co-occurrence networks (SpiecEasi; SparCC option), and functional profiling options

Read data and basic objects

  


#Howmany ESKAPEE pathogens water vs plastic abundance

#Define ESKAPEE list
eskapee_list <- c(
  "Enterococcus faecium",
  "Staphylococcus aureus",
  "Klebsiella pneumoniae",
  "Acinetobacter baumannii",
  "Pseudomonas aeruginosa",
  "Enterobacter",
  "Escherichia coli"
)

#identifyESKAPEE species in the dataset

library(stringr)

eskapee_detected <- rownames(counts_mat)[
  str_detect(rownames(counts_mat),
             paste(eskapee_list, collapse="|"))
]

eskapee_detected
length(eskapee_detected)

## ---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)

## ---- 1) Define ESKAPEE patterns ----
eskapee_list <- c(
  "Enterococcus faecium",
  "Staphylococcus aureus",
  "Klebsiella pneumoniae",
  "Acinetobacter baumannii",
  "Pseudomonas aeruginosa",
  "Enterobacter",
  "Escherichia coli"
)

## ---- 2) Detect ESKAPEE taxa in your species table ----
# counts_mat rows are species names (as in your output)
eskapee_detected <- rownames(counts_mat)[
  str_detect(rownames(counts_mat), paste(eskapee_list, collapse = "|"))
]

# Remove non-clinical "endosymbiont" hits
eskapee_clean <- eskapee_detected[!grepl("endosymbiont", eskapee_detected)]

cat("ESKAPEE detected (clean):", length(eskapee_clean), "\n")
print(eskapee_clean)

## ---- 3) Define Plastic vs Water sample IDs from metadata ----
plastic_samples <- rownames(meta2)[meta2$SampleType == "Plastic"]
water_samples   <- rownames(meta2)[meta2$SampleType == "Water"]

## ---- 4) Per-sample ESKAPEE abundance (RAW) ----
esk_raw_per_sample <- data.frame(
  SampleID = c(plastic_samples, water_samples),
  Habitat  = c(rep("Plastic", length(plastic_samples)),
               rep("Water",   length(water_samples))),
  Eskapee_raw = c(
    colSums(counts_mat[eskapee_clean, plastic_samples, drop = FALSE]),
    colSums(counts_mat[eskapee_clean, water_samples, drop = FALSE])
  )
)

## ---- 5) Summary mean ± SE ----
esk_sum <- esk_raw_per_sample %>%
  group_by(Habitat) %>%
  summarise(
    mean = mean(Eskapee_raw),
    se   = sd(Eskapee_raw) / sqrt(n()),
    .groups = "drop"
  )

print(esk_sum)

## ---- 6) Wilcoxon rank-sum test ----
wilcox_test <- wilcox.test(Eskapee_raw ~ Habitat,
                           data = esk_raw_per_sample,
                           exact = FALSE)
pval <- wilcox_test$p.value
cat("Wilcoxon p-value:", pval, "\n")

## ---- 7) Significance stars ----
stars <- ifelse(pval < 0.001, "***",
         ifelse(pval < 0.01,  "**",
         ifelse(pval < 0.05,  "*", "ns")))

## ---- 8) Plot (bar + SE + bracket + stars) ----
y_max <- max(esk_sum$mean + esk_sum$se)

ggplot(esk_sum, aes(x = Habitat, y = mean, fill = Habitat)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.15, linewidth = 1) +
  scale_fill_manual(values = c("Plastic" = "#8B0000",
                               "Water"   = "#00008B")) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1.3),
    axis.ticks = element_line(color = "black", linewidth = 1.2),
    legend.position = "none"
  ) +
  # Bracket
  geom_segment(aes(x = 1, xend = 2, y = y_max*1.08, yend = y_max*1.08),
               linewidth = 1) +
  geom_segment(aes(x = 1, xend = 1, y = y_max*1.06, yend = y_max*1.08),
               linewidth = 1) +
  geom_segment(aes(x = 2, xend = 2, y = y_max*1.06, yend = y_max*1.08),
               linewidth = 1) +
  # Stars
  annotate("text", x = 1.5, y = y_max*1.11, label = stars, size = 8) +
  labs(
    title = "ESKAPEE abundance per sample",
    y = "Mean total read count (±SE)",
    x = ""
  )


#Network Analysis

## ============================================================
## Plastisphere Co-occurrence Network: Pie nodes by Environment
## SpiecEasi | Plastic vs Water | ESKAPEE numbered
## ============================================================

## ---- 0) Libraries ----
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(SpiecEasi)
library(igraph)

## ---- 1) Load data ----
counts_path <- "/home/meti/Downloads/Lig_wf-metagenomics-counts-species.csv"
meta_path   <- "/home/meti/Downloads/Ligation_all_metadata (1).txt"

counts <- read_csv(counts_path, show_col_types = FALSE)
meta2  <- read.delim(meta_path, header = TRUE, stringsAsFactors = FALSE)

species_col <- colnames(counts)[1]
tax_cols    <- intersect(c("kingdom","phylum","class","order","family","genus","tax"),
                         colnames(counts))

counts_samples_only <- counts %>%
  select(-all_of(tax_cols)) %>%
  column_to_rownames(var = species_col)

counts_mat           <- as.matrix(counts_samples_only)
meta2$SampleID       <- trimws(meta2$SampleID)
colnames(counts_mat) <- trimws(colnames(counts_mat))
meta2                <- meta2[match(colnames(counts_mat), meta2$SampleID), ]

if (!all(meta2$SampleID == colnames(counts_mat))) {
  stop("SampleID mismatch. Check metadata vs count table.")
}

## ---- 2) Remove NA samples ----
meta2      <- meta2[!is.na(meta2$SampleType) & !is.na(meta2$Environment), ]
counts_mat <- counts_mat[, meta2$SampleID]
cat("Samples after NA removal:", ncol(counts_mat), "\n")

## ---- 3) Relative abundance + 70% prevalence filter ----
rel_all              <- sweep(counts_mat, 2, colSums(counts_mat), "/")
rel_all[is.na(rel_all)] <- 0
prevalence           <- rowSums(rel_all > 0)
keep_taxa            <- prevalence >= (0.70 * ncol(rel_all))
rel_filtered         <- rel_all[keep_taxa, ]
cat("Filtered taxa x samples:", dim(rel_filtered), "\n")

## ---- 4) Split Plastic vs Water ----
rel_t       <- t(rel_filtered)
plastic_mat <- rel_t[meta2$SampleType == "Plastic", , drop = FALSE]
water_mat   <- rel_t[meta2$SampleType == "Water",   , drop = FALSE]
plastic_mat <- plastic_mat[, colSums(plastic_mat) > 0, drop = FALSE]
water_mat   <- water_mat[,   colSums(water_mat)   > 0, drop = FALSE]
common_taxa <- intersect(colnames(plastic_mat), colnames(water_mat))
plastic_mat <- plastic_mat[, common_taxa, drop = FALSE]
water_mat   <- water_mat[,   common_taxa, drop = FALSE]
cat("Common taxa:", length(common_taxa), "\n")

## ---- 5) SpiecEasi networks (takes a few minutes) ----
set.seed(1)
se_plastic <- spiec.easi(plastic_mat, method = "mb",
                         lambda.min.ratio = 1e-2, nlambda = 20,
                         pulsar.params = list(rep.num = 20))
se_water   <- spiec.easi(water_mat,   method = "mb",
                         lambda.min.ratio = 1e-2, nlambda = 20,
                         pulsar.params = list(rep.num = 20))

adj_plastic <- getRefit(se_plastic)
adj_water   <- getRefit(se_water)
ig_plastic  <- graph_from_adjacency_matrix(adj_plastic, mode = "undirected", diag = FALSE)
ig_water    <- graph_from_adjacency_matrix(adj_water,   mode = "undirected", diag = FALSE)
V(ig_plastic)$name  <- common_taxa
V(ig_water)$name    <- common_taxa
V(ig_plastic)$label <- NA
V(ig_water)$label   <- NA
cat("Nodes - Plastic:", vcount(ig_plastic), "| Water:", vcount(ig_water), "\n")
cat("Edges - Plastic:", ecount(ig_plastic), "| Water:", ecount(ig_water), "\n")

## ---- 6) Canonical ESKAPEE (strict exact match) ----
eskapee_strict <- c(
  "Enterococcus faecium",
  "Staphylococcus aureus",
  "Klebsiella pneumoniae",
  "Acinetobacter baumannii",
  "Pseudomonas aeruginosa",
  "Enterobacter cloacae",
  "Escherichia coli"
)
eskapee_clean <- common_taxa[common_taxa %in% eskapee_strict]
cat("ESKAPEE nodes found:", length(eskapee_clean), "\n")
print(eskapee_clean)

## ---- 7) Numbered labels (only on ESKAPEE nodes) ----
esk_labels_plastic <- ifelse(V(ig_plastic)$name %in% eskapee_clean,
                              as.character(match(V(ig_plastic)$name, eskapee_clean)), "")
esk_labels_water   <- ifelse(V(ig_water)$name %in% eskapee_clean,
                              as.character(match(V(ig_water)$name, eskapee_clean)), "")

## ---- 8) Environment colors ----
env_levels <- c("River", "Estuary", "Sea", "Control")
env_colors <- c(River   = "#4472C4",   # blue
                Estuary = "#FF69B4",   # pink
                Sea     = "#70AD47",   # green
                Control = "#FFC000")   # yellow

## ---- 9) Per-taxon pie proportions by Environment ----
pie_vals <- lapply(common_taxa, function(taxon) {
  vals <- sapply(env_levels, function(e) {
    idx <- meta2$Environment == e & !is.na(meta2$Environment)
    if (sum(idx) > 0) mean(rel_t[idx, taxon], na.rm = TRUE) else 0
  })
  if (sum(vals) == 0) vals[] <- 1 / length(env_levels)
  vals / sum(vals)
})

## ---- 10) Node and edge styling ----
V(ig_plastic)$size  <- 5
V(ig_water)$size    <- 5
E(ig_plastic)$color <- "grey30"
E(ig_water)$color   <- "grey30"
E(ig_plastic)$width <- 0.6
E(ig_water)$width   <- 0.6

## ---- 11) Layout (same seed for both networks) ----
set.seed(1)
lay <- layout_with_fr(ig_plastic)

## ---- 12) Build single-line caption ----
caption_items <- paste0(seq_along(eskapee_clean), ". ", eskapee_clean)
caption_line  <- paste(caption_items, collapse = "     ")

## ---- 13) Save plot ----
png("/home/meti/Downloads/Plastisphere_Network_Pie.png",
    width = 4800, height = 2400, res = 300)

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(layout_matrix, heights = c(6, 1))

# --- Plastic network ---
par(mar = c(1, 1, 3, 1), bg = "white")
plot(ig_plastic,
     layout             = lay,
     vertex.shape       = "pie",
     vertex.pie         = pie_vals,
     vertex.pie.color   = list(env_colors),
     vertex.size        = V(ig_plastic)$size,
     vertex.label       = esk_labels_plastic,
     vertex.label.cex   = 0.6,
     vertex.label.color = "black",
     vertex.label.font  = 2,
     vertex.frame.color = ifelse(V(ig_plastic)$name %in% eskapee_clean, "#8B0000", "grey40"),
     vertex.frame.width = ifelse(V(ig_plastic)$name %in% eskapee_clean, 2.5, 0.5),
     edge.color         = E(ig_plastic)$color,
     edge.width         = E(ig_plastic)$width,
     main               = "Plastisphere Network (Plastic)")

legend("topleft",
       legend = env_levels,
       fill   = env_colors[env_levels],
       title  = "Environment",
       bty    = "n", cex = 0.8)

# --- Water network ---
par(mar = c(1, 1, 3, 1), bg = "white")
plot(ig_water,
     layout             = lay,
     vertex.shape       = "pie",
     vertex.pie         = pie_vals,
     vertex.pie.color   = list(env_colors),
     vertex.size        = V(ig_water)$size,
     vertex.label       = esk_labels_water,
     vertex.label.cex   = 0.6,
     vertex.label.color = "black",
     vertex.label.font  = 2,
     vertex.frame.color = ifelse(V(ig_water)$name %in% eskapee_clean, "#8B0000", "grey40"),
     vertex.frame.width = ifelse(V(ig_water)$name %in% eskapee_clean, 2.5, 0.5),
     edge.color         = E(ig_water)$color,
     edge.width         = E(ig_water)$width,
     main               = "Water Column Network (Water)")

# --- Caption panel (single line using mtext) ---
par(mar = c(0, 0, 0, 0), bg = "white")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.75, "ESKAPEE Pathogens detected in network:",
     cex = 0.9, font = 2, adj = 0.5)
mtext(caption_line, side = 1, line = -3, cex = 0.65, font = 3, adj = 0.5)

dev.off()
cat("Saved to: /home/meti/Downloads/Plastisphere_Network_Pie.png\n")

## ---- 14) Network statistics ----
network_stats <- data.frame(
  Network         = c("Plastic", "Water"),
  Nodes           = c(vcount(ig_plastic), vcount(ig_water)),
  Edges           = c(ecount(ig_plastic), ecount(ig_water)),
  MeanDegree      = round(c(mean(degree(ig_plastic)), mean(degree(ig_water))), 2),
  ClusteringCoeff = round(c(transitivity(ig_plastic, type = "global"),
                             transitivity(ig_water,   type = "global")), 3),
  Modularity      = round(c(modularity(cluster_fast_greedy(ig_plastic)),
                             modularity(cluster_fast_greedy(ig_water))), 3)
)
print(network_stats)



##(Matan’s comment to show them in different color)
# 1st we need to filter them with this code

# --- Rebuild rel_filtered from scratch ---
counts_path <- "C:/Users/metig/Downloads/Lig_wf-metagenomics-counts-species.csv"
meta_path   <- "C:/Users/metig/Downloads/Ligation_all_metadata (1).txt"

counts_raw <- read_csv(counts_path, show_col_types = FALSE)
meta2      <- read.delim(meta_path, header = TRUE, stringsAsFactors = FALSE)

# Species column
species_col <- colnames(counts_raw)[1]

# Remove taxonomy columns
tax_cols2 <- intersect(c("total","kingdom","phylum","class","order",
                          "family","genus","tax","superkingdom"),
                        colnames(counts_raw))

counts_mat <- counts_raw %>%
  select(-all_of(tax_cols2)) %>%
  column_to_rownames(var = species_col) %>%
  as.matrix()

# Match metadata
meta2$SampleID <- trimws(meta2$SampleID)
colnames(counts_mat) <- trimws(colnames(counts_mat))
meta2 <- meta2[match(colnames(counts_mat), meta2$SampleID), ]

# Relative abundance
rel_all <- sweep(counts_mat, 2, colSums(counts_mat), "/")
rel_all[is.na(rel_all)] <- 0

# Prevalence filter 70%
prevalence <- rowSums(rel_all > 0)
keep_taxa  <- prevalence >= (0.70 * ncol(rel_all))
rel_filtered <- rel_all[keep_taxa, ]
cat("Filtered taxa:", nrow(rel_filtered), "\n")

# Split plastic vs water
rel_t <- t(rel_filtered)
plastic_mat <- rel_t[meta2$SampleType == "Plastic", , drop = FALSE]
water_mat   <- rel_t[meta2$SampleType == "Water",   , drop = FALSE]
plastic_mat <- plastic_mat[, colSums(plastic_mat) > 0, drop = FALSE]
water_mat   <- water_mat[,   colSums(water_mat)   > 0, drop = FALSE]
common_taxa <- intersect(colnames(plastic_mat), colnames(water_mat))
plastic_mat <- plastic_mat[, common_taxa, drop = FALSE]
water_mat   <- water_mat[,   common_taxa, drop = FALSE]
cat("Common taxa:", length(common_taxa), "\n")

# Run SpiecEasi
set.seed(1)
se_plastic <- spiec.easi(plastic_mat,
                          method = "mb",
                          lambda.min.ratio = 1e-2,
                          nlambda = 20,
                          pulsar.params = list(rep.num = 20))

adj_plastic <- getRefit(se_plastic)
ig_plastic  <- graph_from_adjacency_matrix(adj_plastic,
                                            mode = "undirected", diag = FALSE)
V(ig_plastic)$name <- common_taxa

# ESKAPEE nodes
eskapee_list <- c("Enterococcus faecium", "Staphylococcus aureus",
                  "Klebsiella pneumoniae", "Acinetobacter baumannii",
                  "Pseudomonas aeruginosa", "Enterobacter", "Escherichia coli")
eskapee_detected <- common_taxa[str_detect(common_taxa,
                                            paste(eskapee_list, collapse="|"))]
eskapee_clean <- eskapee_detected[!grepl("endosymbiont", eskapee_detected)]
cat("ESKAPEE nodes:", length(eskapee_clean), "\n")
print(eskapee_clean)

cat("All objects ready! Now run the pie network plot.\n")
 
###then build the pie network 

# Run water network too
set.seed(1)
se_water <- spiec.easi(water_mat,
                       method = "mb",
                       lambda.min.ratio = 1e-2,
                       nlambda = 20,
                       pulsar.params = list(rep.num = 20))

adj_water <- getRefit(se_water)
ig_water  <- graph_from_adjacency_matrix(adj_water,
                                          mode = "undirected", diag = FALSE)
V(ig_water)$name <- common_taxa

# Same pie values for water network
pie_values_water <- lapply(V(ig_water)$name, function(nm) {
  as.numeric(env_abund_norm[nm, env_order])
})

deg_water       <- degree(ig_water)
node_size_water <- rep(4, vcount(ig_water))
node_size_water[V(ig_water)$name %in% eskapee_clean] <- 12
node_size_water[deg_water > quantile(deg_water, 0.90)] <- 
  node_size_water[deg_water > quantile(deg_water, 0.90)] + 2

frame_col_water <- ifelse(V(ig_water)$name %in% eskapee_clean, "#8B0000", "#333333")
frame_wid_water <- ifelse(V(ig_water)$name %in% eskapee_clean, 3, 0.3)

E(ig_water)$color <- "#AAAAAA"
E(ig_water)$width <- 0.5

# Same layout for fair comparison
set.seed(42)
lay_plastic <- layout_with_fr(ig_plastic)
lay_water   <- layout_with_fr(ig_water)

# --- Plot side by side ---
windows(width = 22, height = 11)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2), bg = "white")

# Plastic
plot(ig_plastic,
     layout             = lay_plastic,
     vertex.shape       = "pie",
     vertex.pie         = pie_values,
     vertex.pie.color   = list(env_colors[env_order]),
     vertex.size        = node_size,
     vertex.label       = ifelse(V(ig_plastic)$name %in% eskapee_clean,
                                  V(ig_plastic)$name, NA),
     vertex.label.cex   = 0.5,
     vertex.label.color = "black",
     vertex.label.font  = 3,
     vertex.label.dist  = 1.2,
     vertex.frame.color = frame_col,
     vertex.frame.width = frame_wid,
     edge.color         = "#AAAAAA",
     edge.width         = 0.5,
     main               = "Plastisphere Network (Plastic)"
)

# Water
plot(ig_water,
     layout             = lay_water,
     vertex.shape       = "pie",
     vertex.pie         = pie_values_water,
     vertex.pie.color   = list(env_colors[env_order]),
     vertex.size        = node_size_water,
     vertex.label       = ifelse(V(ig_water)$name %in% eskapee_clean,
                                  V(ig_water)$name, NA),
     vertex.label.cex   = 0.5,
     vertex.label.color = "black",
     vertex.label.font  = 3,
     vertex.label.dist  = 1.2,
     vertex.frame.color = frame_col_water,
     vertex.frame.width = frame_wid_water,
     edge.color         = "#AAAAAA",
     edge.width         = 0.5,
     main               = "Water Column Network (Water)"
)

# Shared legend
legend("bottomleft",
       legend = env_order,
       pch    = 21,
       pt.bg  = env_colors[env_order],
       col    = "black",
       pt.cex = 1.8,
       bty    = "n",
       cex    = 0.8,
       title  = "Environment"
)
legend("bottomright",
       legend = c("ESKAPEE pathogen", "Other taxon"),
       pch    = 21,
       pt.bg  = c("#8B0000", "#555555"),
       col    = c("#8B0000", "#333333"),
       pt.cex = c(2.5, 1.5),
       bty    = "n",
       cex    = 0.8,
       title  = "Node type"
)

# Save
dev.copy(png, "D:/network_pie_plastic_vs_water.png",
         width = 6400, height = 3200, res = 300)
dev.off()
cat("Saved!\n")


Shared and unique taxa

# ── Create meta2 from your existing meta object ───────────────────────────────
meta2 <- meta

# Add Environment column based on Site
meta2$Environment <- dplyr::recode(meta2$Site,
  "R2" = "River",
  "R3" = "Estuary",
  "S1" = "Sea",
  "Sc" = "Control"
)

# Verify
cat("SampleType values:", unique(meta2$SampleType), "\n")
cat("Environment values:", unique(meta2$Environment), "\n")
cat("Dimensions:", nrow(meta2), "samples\n")
print(head(meta2))



install.packages("VennDiagram")
library(VennDiagram)
library(grid)
presence_cutoff <- 0   # change to 4 for >=5 reads

get_taxa_set <- function(sample_ids) {
  rownames(counts_mat)[rowSums(counts_mat[, sample_ids, drop=FALSE] > presence_cutoff) > 0]
}

# ---- WATER sets (by location) ----
water_meta <- meta2[meta2$SampleType == "Water", ]

water_sets <- list(
  River   = get_taxa_set(water_meta$SampleID[water_meta$Environment == "River"]),
  Estuary = get_taxa_set(water_meta$SampleID[water_meta$Environment == "Estuary"]),
  Sea     = get_taxa_set(water_meta$SampleID[water_meta$Environment == "Sea"]),
  Control = get_taxa_set(water_meta$SampleID[water_meta$Environment == "Control"])
)

# ---- PLASTIC sets (by location) ----
plastic_meta <- meta2[meta2$SampleType == "Plastic", ]

plastic_sets <- list(
  River   = get_taxa_set(plastic_meta$SampleID[plastic_meta$Environment == "River"]),
  Estuary = get_taxa_set(plastic_meta$SampleID[plastic_meta$Environment == "Estuary"]),
  Sea     = get_taxa_set(plastic_meta$SampleID[plastic_meta$Environment == "Sea"]),
  Control = get_taxa_set(plastic_meta$SampleID[plastic_meta$Environment == "Control"])
)

# Quick check sizes
sapply(water_sets, length)
sapply(plastic_sets, length)

venn_water <- venn.diagram(
  x = water_sets,
  filename = NULL,
  fill = c("#8B0000", "#00008B", "#006400", "#4B0082"), # dark red, dark blue, dark green, indigo
  alpha = 0.55,
  lwd = 2,
  col = "black",
  cex = 1.2,
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.08,
  main = "WATER: shared and unique taxa by location",
  main.cex = 1.4
)

venn_plastic <- venn.diagram(
  x = plastic_sets,
  filename = NULL,
  fill = c("#8B0000", "#00008B", "#006400", "#4B0082"),
  alpha = 0.55,
  lwd = 2,
  col = "black",
  cex = 1.2,
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.08,
  main = "PLASTIC: shared and unique taxa by location",
  main.cex = 1.4
)

# Side-by-side plotting
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(venn_water)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(venn_plastic)
popViewport(2)

# PNG
png("Venn_Water_vs_Plastic_byLocation.png", width=2600, height=1400, res=300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(venn_water)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(venn_plastic)
popViewport(2)
dev.off()

# PDF
pdf("Venn_Water_vs_Plastic_byLocation.pdf", width=12, height=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(venn_water)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(venn_plastic)
popViewport(2)
dev.off()

SOURCE TRACKING

# ============================================================
# FEAST — Taninim River Plastisphere
# Rebuilt from working March 9 session
# Run in RStudio on Windows
# ============================================================

# --- Step 1: Install packages ---
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("cozygene/FEAST", dependencies = TRUE, quiet = TRUE)

for (pkg in c("ggplot2", "dplyr", "tidyr", "readr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(FEAST)
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
cat("Libraries loaded\n")

# --- Step 2: Load data ---
counts_raw <- read_csv("C:/Users/metig/Downloads/Lig_wf-metagenomics-counts-species.csv",
                       show_col_types = FALSE)
meta       <- read.delim("C:/Users/metig/Downloads/Ligation_all_metadata (1) (1).txt",
                         header = TRUE, stringsAsFactors = FALSE)

# --- Step 3: Fix column names ---
# Handle #SampleID / X.SampleID
colnames(meta) <- gsub("^#", "", colnames(meta))
colnames(meta) <- gsub("^X\\.", "", colnames(meta))
meta$SampleID  <- trimws(meta$SampleID)
cat(sprintf("Metadata columns: %s\n", paste(colnames(meta), collapse = ", ")))

# --- Step 4: Build count matrix ---
tax_cols <- intersect(
  c("total","kingdom","phylum","class","order","family","genus","tax","superkingdom"),
  colnames(counts_raw)
)

count_mat_raw <- counts_raw %>%
  select(-all_of(tax_cols)) %>%
  column_to_rownames(var = "species") %>%
  as.matrix()

# FIX: S1P1 -> S1P4
if ("S1P1" %in% colnames(count_mat_raw)) {
  colnames(count_mat_raw)[colnames(count_mat_raw) == "S1P1"] <- "S1P4"
  cat("Fixed: S1P1 -> S1P4\n")
}

# Transpose: samples x species
count_mat <- t(count_mat_raw)
cat(sprintf("count_mat dimensions: %d samples x %d species\n",
            nrow(count_mat), ncol(count_mat)))
cat("Row names:\n"); print(rownames(count_mat))

# --- Step 5: Filter species ---
# Keep species with >=50 total reads AND present in >=2 samples in at least one site
sites <- list(
  R2 = c("R2P1","R2P2","R2P3","R2W1","R2W2","R2W3"),
  R3 = c("R3P1","R3P2","R3P3","R3W1","R3W2","R3W3"),
  S1 = c("S1P4","S1P2","S1P3","S1W1","S1W2","S1W3"),
  Sc = c("ScP1","ScP2","ScP3","ScW1","ScW2","ScW3")
)

pass_total <- colSums(count_mat) >= 50

pass_site <- apply(count_mat, 2, function(sp_counts) {
  any(sapply(sites, function(s) {
    valid <- s[s %in% rownames(count_mat)]
    sum(sp_counts[valid] > 0) >= 2
  }))
})

count_mat <- count_mat[, pass_total & pass_site]
cat(sprintf("After filtering: %d species retained\n", ncol(count_mat)))

# --- Step 6: Subsampling function (replaces COVERAGE parameter) ---
subsample_counts <- function(mat, depth = 10000) {
  t(apply(mat, 1, function(x) {
    if (sum(x) <= depth) return(x)
    as.integer(rmultinom(1, depth, x))
  }))
}

# --- Step 7: FEAST runner function ---
run_feast <- function(meta_feast, label, outdir) {
  samps <- rownames(meta_feast)

  # Check all samples exist
  missing <- samps[!samps %in% rownames(count_mat)]
  if (length(missing) > 0) {
    cat(sprintf("WARNING: missing samples: %s\n", paste(missing, collapse = ", ")))
    return(NULL)
  }

  mat <- count_mat[samps, ]
  mat[is.na(mat)] <- 0
  mat <- mat[, colSums(mat) > 0]
  mat <- round(mat)
  storage.mode(mat) <- "integer"

  set.seed(42)
  mat <- subsample_counts(mat, depth = 10000)
  colnames(mat) <- colnames(count_mat[samps, colSums(count_mat[samps,]) > 0])

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("Running %s | %d samples x %d OTUs\n", label, nrow(mat), ncol(mat)))

  FEAST(C                      = mat,
        metadata               = meta_feast,
        different_sources_flag = 0,
        dir_path               = outdir,
        outfile                = label)
}

# --- Step 8: Define metadata for all 4 scenarios ---

# Scenario 1: R3P Estuary Plastic — sources: R2P + R2W + R3W
feast_meta_sc1 <- data.frame(
  Env        = c(rep("R2 Plastisphere",3), rep("R2 Water",3), rep("R3 Water",3),
                 rep("R3 Plastisphere",3)),
  SourceSink = c(rep("Source",9), rep("Sink",3)),
  id         = c(1,1,1, 1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R2W1","R2W2","R2W3",
                 "R3W1","R3W2","R3W3",
                 "R3P1","R3P2","R3P3"),
  stringsAsFactors = FALSE
)

# Scenario 2: R3W Estuary Water — sources: R2P + R2W
feast_meta_sc2 <- data.frame(
  Env        = c(rep("R2 Plastisphere",3), rep("R2 Water",3),
                 rep("R3 Water",3)),
  SourceSink = c(rep("Source",6), rep("Sink",3)),
  id         = c(1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R2W1","R2W2","R2W3",
                 "R3W1","R3W2","R3W3"),
  stringsAsFactors = FALSE
)

# Scenario 3: S1P Sea Plastic — sources: R2P + R2W + R3P + R3W + S1W
feast_meta_sc3 <- data.frame(
  Env        = c(rep("R2 Plastisphere",3), rep("R2 Water",3),
                 rep("R3 Plastisphere",3), rep("R3 Water",3),
                 rep("S1 Water",3),
                 rep("S1 Plastisphere",3)),
  SourceSink = c(rep("Source",15), rep("Sink",3)),
  id         = c(1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R2W1","R2W2","R2W3",
                 "R3P1","R3P2","R3P3",
                 "R3W1","R3W2","R3W3",
                 "S1W1","S1W2","S1W3",
                 "S1P4","S1P2","S1P3"),
  stringsAsFactors = FALSE
)

# Scenario 4: S1W Sea Water — sources: R2P + R3P + R3W
feast_meta_sc4 <- data.frame(
  Env        = c(rep("R2 Plastisphere",3), rep("R3 Plastisphere",3),
                 rep("R3 Water",3),
                 rep("S1 Water",3)),
  SourceSink = c(rep("Source",9), rep("Sink",3)),
  id         = c(1,1,1, 1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R3P1","R3P2","R3P3",
                 "R3W1","R3W2","R3W3",
                 "S1W1","S1W2","S1W3"),
  stringsAsFactors = FALSE
)

# --- Step 9: Run all 4 scenarios ---
sc1 <- run_feast(feast_meta_sc1, "sc1", "D:/feast_sc1/")
sc2 <- run_feast(feast_meta_sc2, "sc2", "D:/feast_sc2/")
sc3 <- run_feast(feast_meta_sc3, "sc3", "D:/feast_sc3/")
sc4 <- run_feast(feast_meta_sc4, "sc4", "D:/feast_sc4/")

cat("\n=== Scenario 1: R3P Estuary Plastic ===\n")
print(round(sc1$data_prop * 100, 2))

cat("\n=== Scenario 2: R3W Estuary Water ===\n")
print(round(sc2$data_prop * 100, 2))

cat("\n=== Scenario 3: S1P Sea Plastic ===\n")
print(round(sc3$data_prop * 100, 2))

cat("\n=== Scenario 4: S1W Sea Water ===\n")
print(round(sc4$data_prop * 100, 2))

# --- Step 10: Summarize and plot ---
summarize_scenario <- function(data_prop, sink_label) {
  df <- as.data.frame(data_prop[, "pred_emnoise_all", drop = FALSE])
  df$source_replicate <- rownames(df)
  df$Env <- gsub("_[^_]+$", "", df$source_replicate)
  env_sums <- df %>%
    group_by(Env) %>%
    summarise(Proportion = sum(pred_emnoise_all) * 100, .groups = "drop")
  env_sums$Proportion <- env_sums$Proportion / sum(env_sums$Proportion) * 100
  env_sums$Sink <- sink_label
  return(env_sums)
}

sc1_sum <- summarize_scenario(sc1$data_prop, "R3 Plastisphere (Estuary)")
sc2_sum <- summarize_scenario(sc2$data_prop, "R3 Water (Estuary)")
sc3_sum <- summarize_scenario(sc3$data_prop, "S1 Plastisphere (Sea)")
sc4_sum <- summarize_scenario(sc4$data_prop, "S1 Water (Sea)")

plot_df <- bind_rows(sc1_sum, sc2_sum, sc3_sum, sc4_sum)

plot_df$Sink <- factor(plot_df$Sink, levels = c(
  "R3 Plastisphere (Estuary)",
  "R3 Water (Estuary)",
  "S1 Plastisphere (Sea)",
  "S1 Water (Sea)"
))

source_order <- c("R2 Plastisphere","R3 Plastisphere",
                  "R2 Water","R3 Water","S1 Water","unknown")
plot_df$Env <- factor(plot_df$Env, levels = source_order)

plot_df$Label <- ifelse(plot_df$Proportion >= 3,
                        paste0(round(plot_df$Proportion, 1), "%"), "")

source_colors <- c(
  "R2 Plastisphere" = "#7B2D8B",
  "R3 Plastisphere" = "#2B5EA7",
  "R2 Water"        = "#E8638C",
  "R3 Water"        = "#6B1A1A",
  "S1 Water"        = "#148F77",
  "unknown"         = "#D4860A"
)

p_st <- ggplot(plot_df, aes(x = Sink, y = Proportion, fill = Env)) +
  geom_bar(stat = "identity", width = 0.55,
           color = "black", linewidth = 0.8) +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = source_colors, name = "Source") +
  scale_x_discrete(limits = rev(c(
    "R3 Plastisphere (Estuary)",
    "R3 Water (Estuary)",
    "S1 Plastisphere (Sea)",
    "S1 Water (Sea)"
  ))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101),
                     breaks = seq(0, 100, 20),
                     labels = paste0(seq(0, 100, 20), "%")) +
  coord_flip() +
  labs(
    title    = "Microbial Source Tracking of Taninim Plastisphere Communities",
    subtitle = "FEAST Bayesian source tracking | Taninim River-Mediterranean System",
    x        = NULL,
    y        = "Source Contribution (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 10, color = "grey40"),
    axis.text.y       = element_text(size = 11, color = "black", face = "bold"),
    axis.text.x       = element_text(size = 11, color = "black", face = "bold"),
    axis.title.x      = element_text(size = 12, color = "black", face = "bold"),
    axis.line         = element_line(color = "black", linewidth = 1.2),
    axis.ticks        = element_line(color = "black", linewidth = 1.0),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title      = element_text(face = "bold", size = 11),
    legend.text       = element_text(size = 10),
    legend.position   = "right",
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

print(p_st)
ggsave("D:/FEAST_sourcetracking_4scenarios.png", p_st,
       width = 10, height = 6, dpi = 300, bg = "white")

write.csv(plot_df, "D:/FEAST_summary.csv", row.names = FALSE)
cat("Saved: D:/FEAST_sourcetracking_4scenarios.png\n")
cat("Saved: D:/FEAST_summary.csv\n")
cat("\nDone!\n")



Python tracking the source

# ============================================================
# SourceTracker2 Analysis: Taninim River Plastisphere
# FIXED VERSION - proper Gibbs sampling convergence
# Copy-paste ALL of this into ONE Colab cell and run
# ============================================================

# --- Step 1: Upload files ---
from google.colab import files
print("Upload your 2 files now:")
uploaded = files.upload()

# --- Step 2: Install and import ---
!pip install scikit-bio -q

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# --- Step 3: Load data ---
counts_df = pd.read_csv('Lig_wf-metagenomics-counts-species.csv')
meta = pd.read_csv('Ligation_all_metadata (1).txt', sep='\t')

counts_df.columns = counts_df.columns.str.strip().str.replace('\r', '')
meta.columns = meta.columns.str.strip().str.replace('\r', '')
meta['SampleID'] = meta['SampleID'].str.strip().str.replace('\r', '')

# Rename S1P1 to S1P4 in count table to match metadata
counts_df.columns = counts_df.columns.str.replace('S1P1', 'S1P4')

tax_cols = ['species','kingdom','phylum','class','order','family','genus','tax']
sample_cols = [c for c in counts_df.columns if c not in tax_cols]
print(f'Detected {len(sample_cols)} sample columns: {sample_cols}')

count_matrix = counts_df.set_index('species')[sample_cols].T
meta = meta.set_index('SampleID').loc[sample_cols]

print(f'Count matrix: {count_matrix.shape[0]} samples x {count_matrix.shape[1]} species')
display(meta)

# --- Step 4: Bayesian Source Tracker (FIXED) ---
def gibbs_source_tracker(source_envs, sink_samples, count_matrix,
                          alpha1=0.001, alpha2=0.01, n_iter=50, n_burn=25,
                          rarefaction_depth=91000):
    """
    Bayesian source tracking via Gibbs sampling (Knights et al. 2011).
    Vectorized implementation for large rarefaction depths.
    """
    source_names = list(source_envs.keys()) + ['Unknown']
    n_sources = len(source_names)

    all_samples_list = []
    for samples in source_envs.values():
        all_samples_list.extend(samples)
    all_samples_list.extend(sink_samples)

    cm = count_matrix.loc[all_samples_list].copy()
    species = cm.columns[cm.sum() > 0]
    cm = cm[species]
    n_species = len(species)

    # Source probability distributions (with Dirichlet prior)
    source_dists = np.zeros((n_sources, n_species))
    for i, (env_name, samples) in enumerate(source_envs.items()):
        source_counts = cm.loc[samples].sum(axis=0).values + alpha2
        source_dists[i] = source_counts / source_counts.sum()
    source_dists[-1] = 1.0 / n_species

    # Precompute log source distributions for numerical stability
    log_source_dists = np.log(source_dists + 1e-300)

    results = {}

    for sink_id in sink_samples:
        print(f'  Processing {sink_id}...')
        sink_counts = cm.loc[sink_id].values.astype(float)
        total_reads = sink_counts.sum()

        if total_reads == 0:
            results[sink_id] = {name: 0.0 for name in source_names}
            continue

        # Rarefy sink to rarefaction_depth
        sink_probs = sink_counts / total_reads
        rng = np.random.default_rng(42)
        depth = min(rarefaction_depth, int(total_reads))
        read_species = rng.choice(n_species, size=depth, p=sink_probs)
        n_reads = len(read_species)

        # Precompute source likelihoods for each read: shape (n_reads, n_sources)
        # source_dists[:, read_species[r]] gives likelihood of read r under each source
        read_likelihoods = source_dists[:, read_species].T  # shape: (n_reads, n_sources)

        # Initialize assignments using source likelihoods
        init_probs = read_likelihoods / read_likelihoods.sum(axis=1, keepdims=True)
        assignments = np.array([rng.choice(n_sources, p=init_probs[r]) for r in range(n_reads)])
        source_counts_tracker = np.bincount(assignments, minlength=n_sources).astype(float)

        collected = []

        for iteration in range(n_iter):
            # Shuffle read order for better mixing
            order = rng.permutation(n_reads)

            for r in order:
                old_source = assignments[r]
                source_counts_tracker[old_source] -= 1

                # Vectorized conditional probability
                probs = (source_counts_tracker + alpha1) * read_likelihoods[r]
                probs /= probs.sum()

                new_source = rng.choice(n_sources, p=probs)
                assignments[r] = new_source
                source_counts_tracker[new_source] += 1

            if iteration >= n_burn:
                proportions = (source_counts_tracker + alpha1) / (n_reads + alpha1 * n_sources)
                collected.append(proportions.copy())

            if (iteration + 1) % 5 == 0:
                current = source_counts_tracker / n_reads * 100
                status = ', '.join([f'{source_names[s]}: {current[s]:.1f}%' for s in range(n_sources)])
                print(f'    Iter {iteration+1}/{n_iter}: {status}')

        mean_props = np.mean(collected, axis=0)
        mean_props /= mean_props.sum()
        results[sink_id] = {name: prop * 100 for name, prop in zip(source_names, mean_props)}
        print(f'  Done: {", ".join([f"{k}: {v:.1f}%" for k,v in results[sink_id].items()])}')

    return pd.DataFrame(results).T, source_names

print('Source tracker ready (fixed version).')

# --- Step 5: Run 3 scenarios ---
np.random.seed(42)

print('=' * 60)
print('Scenario 1: R3 Plastisphere (Estuary) origins')
print('=' * 60)
sc1_results, sc1_sources = gibbs_source_tracker(
    source_envs={
        'R2 Plastisphere': ['R2P1','R2P2','R2P3'],
        'R2 Water': ['R2W1','R2W2','R2W3'],
        'R3 Water': ['R3W1','R3W2','R3W3'],
    },
    sink_samples=['R3P1','R3P2','R3P3'],
    count_matrix=count_matrix
)
print('\n--- Scenario 1 Results ---')
print(sc1_results.round(2))
print(f'\nMean:\n{sc1_results.mean().round(2)}')

print('\n' + '=' * 60)
print('Scenario 2: S1 Plastisphere (Sea) origins')
print('=' * 60)
sc2_results, sc2_sources = gibbs_source_tracker(
    source_envs={
        'R2 Plastisphere': ['R2P1','R2P2','R2P3'],
        'R3 Plastisphere': ['R3P1','R3P2','R3P3'],
        'S1 Water': ['S1W1','S1W2','S1W3'],
        'R3 Water': ['R3W1','R3W2','R3W3'],
    },
    sink_samples=['S1P4','S1P2','S1P3'],
    count_matrix=count_matrix
)
print('\n--- Scenario 2 Results ---')
print(sc2_results.round(2))
print(f'\nMean:\n{sc2_results.mean().round(2)}')

print('\n' + '=' * 60)
print('Scenario 3: S1 Water (Sea) origins')
print('=' * 60)
sc3_results, sc3_sources = gibbs_source_tracker(
    source_envs={
        'R2 Water': ['R2W1','R2W2','R2W3'],
        'R3 Water': ['R3W1','R3W2','R3W3'],
        'R2 Plastisphere': ['R2P1','R2P2','R2P3'],
    },
    sink_samples=['S1W1','S1W2','S1W3'],
    count_matrix=count_matrix
)
print('\n--- Scenario 3 Results ---')
print(sc3_results.round(2))
print(f'\nMean:\n{sc3_results.mean().round(2)}')

# --- Step 6: Pie chart figure ---
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

source_colors = {
    'R2 Plastisphere': '#E74C3C',
    'R2 Water': '#3498DB',
    'R3 Plastisphere': '#E67E22',
    'R3 Water': '#2980B9',
    'S1 Water': '#1ABC9C',
    'Unknown': '#7F8C8D'
}

scenarios = [
    ('R3 Plastisphere\n(Estuary)', sc1_results),
    ('S1 Plastisphere\n(Sea)', sc2_results),
    ('S1 Water\n(Sea)', sc3_results),
]

for idx, (title, results_df) in enumerate(scenarios):
    ax = axes[idx]
    mean_props = results_df.mean().sort_values(ascending=False)

    labels = mean_props.index.tolist()
    values = mean_props.values
    colors = [source_colors.get(l, '#BDC3C7') for l in labels]

    wedges, texts, autotexts = ax.pie(
        values, labels=None, colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p > 2 else '',
        pctdistance=0.75, startangle=90,
        textprops={'fontsize': 10}
    )

    legend_labels = [f'{l}: {v:.1f}%' for l, v in zip(labels, values) if v > 0.5]
    legend_wedges = [w for w, v in zip(wedges, values) if v > 0.5]
    ax.legend(legend_wedges, legend_labels, loc='center left',
             bbox_to_anchor=(-0.35, 0.5), fontsize=8, frameon=False)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=15)

plt.suptitle('Source Tracking: Origin of Downstream Microbial Communities\n'
             'Bayesian Gibbs Sampling (Knights et al. 2011)',
             fontsize=14, y=1.05)
plt.tight_layout()
plt.savefig('sourcetracker_results.png', dpi=300, bbox_inches='tight')
plt.show()
print('Pie chart saved.')

# --- Step 7: Stacked bar figure ---
fig, ax = plt.subplots(figsize=(10, 6))

all_sources = set()
for _, df in scenarios:
    all_sources.update(df.columns)
all_sources = sorted(all_sources)

sink_labels = ['R3 Plastisphere\n(Estuary)', 'S1 Plastisphere\n(Sea)', 'S1 Water\n(Sea)']
data_for_bar = []
for _, df in scenarios:
    mean_vals = df.mean()
    row = {src: mean_vals.get(src, 0) for src in all_sources}
    data_for_bar.append(row)

bar_df = pd.DataFrame(data_for_bar, index=sink_labels)
bar_colors = [source_colors.get(s, '#BDC3C7') for s in bar_df.columns]
bar_df.plot(kind='barh', stacked=True, color=bar_colors, ax=ax, edgecolor='white', linewidth=0.5)

ax.set_xlabel('Contribution (%)', fontsize=12)
ax.set_title('Source Tracking: Contribution of Upstream Sources\nto Downstream Communities', fontsize=13)
ax.legend(title='Source', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
ax.set_xlim(0, 100)
plt.tight_layout()
plt.savefig('sourcetracker_barplot.png', dpi=300, bbox_inches='tight')
plt.show()
print('Stacked bar saved.')

# --- Step 8: Export ---
summary = pd.DataFrame({
    'R3 Plastisphere (Estuary)': sc1_results.mean(),
    'S1 Plastisphere (Sea)': sc2_results.mean(),
    'S1 Water (Sea)': sc3_results.mean()
}).round(2)

print('\n' + '=' * 60)
print('SUMMARY: Mean source contributions (%)')
print('=' * 60)
display(summary)
summary.to_csv('sourcetracker_summary.csv')

files.download('sourcetracker_results.png')
files.download('sourcetracker_barplot.png')
files.download('sourcetracker_summary.csv')
print('\nDone! All files downloaded.')

Phylogenetics tree After Matan’s comment

# ══════════════════════════════════════════════════════════════════════════════
# CELL 1 — Mount Drive, install tools, set paths
# ══════════════════════════════════════════════════════════════════════════════
from google.colab import drive
drive.mount('/content/drive', force_remount=True)

import os, subprocess

!apt-get update -qq
!apt-get install -y -q mafft trimal iqtree
!pip install biopython -q

# Verify
for tool in ["mafft", "trimal", "iqtree2", "iqtree"]:
    r = subprocess.run(["which", tool], capture_output=True, text=True)
    print(f"{'✓' if r.returncode==0 else '✗'} {tool}: {r.stdout.strip()}")

try:
    from Bio import Entrez, SeqIO, Phylo
    print("✓ biopython")
except:
    print("✗ biopython")

# Output paths — all saved to Drive
drive_output = "/content/drive/MyDrive/HOTSPOT_data/phylo_tree"
os.makedirs(drive_output, exist_ok=True)

PATHS = {
    "raw_fasta"     : f"{drive_output}/sequences_raw.fasta",
    "aligned_fasta" : f"{drive_output}/sequences_aligned.fasta",
    "trimmed_fasta" : f"{drive_output}/sequences_trimmed.fasta",
    "tree"          : f"{drive_output}/taninim_silva_tree",
    "itol_phylum"   : f"{drive_output}/iTOL_phylum_colors.txt",
    "itol_substrate": f"{drive_output}/iTOL_substrate.txt",
    "itol_eskapee"  : f"{drive_output}/iTOL_ESKAPEE_highlight.txt",
    "failed"        : f"{drive_output}/failed_species.txt"
}
print("\nAll paths set. Ready for Cell 2.")

# ══════════════════════════════════════════════════════════════════════════════
# CELL 2 — Upload species list
# ══════════════════════════════════════════════════════════════════════════════
from google.colab import files

print("Upload species_for_tree_top500.txt")
uploaded = files.upload()

fname = [f for f in uploaded.keys() if "species" in f.lower()][0]
with open(fname) as f:
    species_list = [line.strip().strip('"') for line in f if line.strip()]

species_list = [s for s in species_list if not any(
    x in s.lower() for x in
    ["unknown"," sp.","str.","cf.","aff.",
     "bacterium","endosymbiont","candidatus"]
)]

print(f"✓ Total species: {len(species_list)}")
print("First 10:", species_list[:10])

# ══════════════════════════════════════════════════════════════════════════════
# CELL 3 — Download 16S sequences from NCBI
# ══════════════════════════════════════════════════════════════════════════════
import time
Entrez.email = "meti.gizaw@ariel.ac.il"

def fetch_16S(species_name, retries=3):
    for attempt in range(retries):
        try:
            # Species level
            q = (f'"{species_name}"[Organism] AND 16S[Title] '
                 f'AND rRNA[Title] AND 1000:1800[SLEN]')
            h = Entrez.esearch(db="nucleotide", term=q,
                               retmax=1, sort="relevance")
            r = Entrez.read(h); h.close()
            if r["IdList"]:
                h = Entrez.efetch(db="nucleotide", id=r["IdList"][0],
                                  rettype="fasta", retmode="text")
                seq = SeqIO.read(h, "fasta"); h.close()
                time.sleep(0.35)
                return seq
            # Genus fallback
            genus = species_name.split()[0]
            q2 = (f'"{genus}"[Organism] AND 16S[Title] AND '
                  f'rRNA[Title] AND 1000:1800[SLEN] AND "type strain"[Title]')
            h = Entrez.esearch(db="nucleotide", term=q2,
                               retmax=1, sort="relevance")
            r2 = Entrez.read(h); h.close()
            if r2["IdList"]:
                h = Entrez.efetch(db="nucleotide", id=r2["IdList"][0],
                                  rettype="fasta", retmode="text")
                seq = SeqIO.read(h, "fasta"); h.close()
                time.sleep(0.35)
                return seq
            return None
        except Exception as e:
            print(f"  Attempt {attempt+1} failed: {e}")
            time.sleep(2)
    return None

sequences, found, failed = [], [], []
print(f"Fetching 16S for {len(species_list)} species (~20-30 min)...\n")

for i, species in enumerate(species_list):
    print(f"[{i+1}/{len(species_list)}] {species}", end=" ... ")
    seq = fetch_16S(species)
    if seq:
        clean = species.replace(" ","_").replace("/","_")
        seq.id = seq.name = clean
        seq.description = ""
        sequences.append(seq)
        found.append(species)
        print("✓")
    else:
        failed.append(species)
        print("✗")
    if (i+1) % 50 == 0:
        SeqIO.write(sequences, PATHS["raw_fasta"], "fasta")
        print(f"\n── Progress saved: {len(sequences)} sequences ──\n")

SeqIO.write(sequences, PATHS["raw_fasta"], "fasta")
with open(PATHS["failed"], "w") as f:
    f.write("\n".join(failed))

print(f"\n✓ Found: {len(found)} | ✗ Failed: {len(failed)}")
print(f"Saved: {PATHS['raw_fasta']}")

# ══════════════════════════════════════════════════════════════════════════════
# CELL 4 — Align with MAFFT
# ══════════════════════════════════════════════════════════════════════════════
print("Running MAFFT (~5-10 min)...")

result = subprocess.run(
    ["mafft","--auto","--thread","4",
     "--reorder", PATHS["raw_fasta"]],
    capture_output=True, text=True
)
with open(PATHS["aligned_fasta"], "w") as f:
    f.write(result.stdout)

print(f"✓ Aligned {result.stdout.count('>')} sequences")
print(f"Saved: {PATHS['aligned_fasta']}")

# ══════════════════════════════════════════════════════════════════════════════
# CELL 5 — Trim with trimAl
# ══════════════════════════════════════════════════════════════════════════════
result = subprocess.run(
    ["trimal",
     "-in",  PATHS["aligned_fasta"],
     "-out", PATHS["trimmed_fasta"],
     "-automated1", "-fasta"],
    capture_output=True, text=True
)

if os.path.exists(PATHS["trimmed_fasta"]):
    with open(PATHS["trimmed_fasta"]) as f:
        print(f"✓ Trimmed: {f.read().count('>')} sequences")
else:
    # Skip trimming if trimal failed
    import shutil
    shutil.copy(PATHS["aligned_fasta"], PATHS["trimmed_fasta"])
    print("trimAl not available — using aligned fasta directly")

print(f"Saved: {PATHS['trimmed_fasta']}")

# ══════════════════════════════════════════════════════════════════════════════
# CELL 6 — Build tree with IQ-TREE2
# ══════════════════════════════════════════════════════════════════════════════
print("Building ML tree (~30-60 min)...")

# Try iqtree2 first, fall back to iqtree
iqtree_cmd = "iqtree2" if subprocess.run(
    ["which","iqtree2"], capture_output=True).returncode == 0 else "iqtree"

result = subprocess.run(
    [iqtree_cmd,
     "-s", PATHS["trimmed_fasta"],
     "-m", "GTR+G",
     "-B", "1000",
     "-T", "AUTO",
     "--prefix", PATHS["tree"],
     "-redo"],
    capture_output=True, text=True
)

tree_file = PATHS["tree"] + ".treefile"
if os.path.exists(tree_file):
    print(f"✓ Tree saved: {tree_file}")
else:
    print("✗ Tree not found — error:")
    print(result.stderr[-500:])

# ══════════════════════════════════════════════════════════════════════════════
# CELL 7 — Generate all iTOL annotation files
# ══════════════════════════════════════════════════════════════════════════════
import pandas as pd

print("Upload counts_top500.csv")
uploaded2 = files.upload()
fname2 = [f for f in uploaded2.keys() if "counts" in f.lower()][0]
counts = pd.read_csv(fname2, index_col=0)

plastic_cols = [c for c in counts.columns if "P" in c]
water_cols   = [c for c in counts.columns if "W" in c]

# ── Phylum map ────────────────────────────────────────────────────────────────
phylum_map = {
    "Escherichia":"Proteobacteria","Klebsiella":"Proteobacteria",
    "Pseudomonas":"Proteobacteria","Acinetobacter":"Proteobacteria",
    "Vibrio":"Proteobacteria","Alcanivorax":"Proteobacteria",
    "Salmonella":"Proteobacteria","Burkholderia":"Proteobacteria",
    "Hydrogenophaga":"Proteobacteria","Comamonas":"Proteobacteria",
    "Celeribacter":"Proteobacteria","Aquirufa":"Proteobacteria",
    "Arcobacter":"Proteobacteria","Cereibacter":"Proteobacteria",
    "Paraburkholderia":"Proteobacteria","Isoalcanivorax":"Proteobacteria",
    "Marivivens":"Proteobacteria","Tabrizicola":"Proteobacteria",
    "Abyssibius":"Proteobacteria",
    "Polaribacter":"Bacteroidota","Flavobacterium":"Bacteroidota",
    "Bacteroides":"Bacteroidota",
    "Bacillus":"Firmicutes","Staphylococcus":"Firmicutes",
    "Corynebacterium":"Actinobacteriota","Mycobacterium":"Actinobacteriota",
    "Planktothrix":"Cyanobacteria","Synechococcus":"Cyanobacteria",
    "Prochlorococcus":"Cyanobacteria","Planctomyces":"Planctomycetota"
}

phylum_colors = {
    "Proteobacteria"  :"#E41A1C",
    "Bacteroidota"    :"#377EB8",
    "Firmicutes"      :"#4DAF4A",
    "Actinobacteriota":"#984EA3",
    "Cyanobacteria"   :"#FF7F00",
    "Planctomycetota" :"#A65628",
    "Other"           :"#999999"
}

# ── A: Phylum color strip ─────────────────────────────────────────────────────
lines = ["DATASET_COLORSTRIP","SEPARATOR TAB",
         "DATASET_LABEL\tPhylum","COLOR\t#000000",
         "LEGEND_TITLE\tPhylum",
         "LEGEND_SHAPES\t"  +"\t".join(["1"]*len(phylum_colors)),
         "LEGEND_COLORS\t"  +"\t".join(phylum_colors.values()),
         "LEGEND_LABELS\t"  +"\t".join(phylum_colors.keys()),
         "DATA"]
for sp in counts.index:
    clean  = sp.replace(" ","_").replace("/","_")
    phylum = phylum_map.get(sp.split()[0], "Other")
    lines.append(f"{clean}\t{phylum_colors[phylum]}\t{phylum}")
with open(PATHS["itol_phylum"],"w") as f: f.write("\n".join(lines))
print(f"✓ {PATHS['itol_phylum']}")

# ── B: Substrate color strip ──────────────────────────────────────────────────
lines = ["DATASET_COLORSTRIP","SEPARATOR TAB",
         "DATASET_LABEL\tSubstrate","COLOR\t#000000",
         "LEGEND_TITLE\tSubstrate type",
         "LEGEND_SHAPES\t1\t1",
         "LEGEND_COLORS\t#8B0000\t#00008B",
         "LEGEND_LABELS\tPlastic-bound\tPlanktonic",
         "DATA"]
for sp in counts.index:
    clean = sp.replace(" ","_").replace("/","_")
    pt = counts.loc[sp, plastic_cols].sum()
    wt = counts.loc[sp, water_cols].sum()
    color = "#8B0000" if pt >= wt else "#00008B"
    label = "Plastic-bound" if pt >= wt else "Planktonic"
    lines.append(f"{clean}\t{color}\t{label}")
with open(PATHS["itol_substrate"],"w") as f: f.write("\n".join(lines))
print(f"✓ {PATHS['itol_substrate']}")

# ── C: ESKAPEE highlight ──────────────────────────────────────────────────────
eskapee_genera = ["Enterococcus","Staphylococcus","Klebsiella",
                  "Acinetobacter","Pseudomonas","Enterobacter","Escherichia"]
tree = Phylo.read(PATHS["tree"]+".treefile","newick")
eskapee_tips = [c.name for c in tree.find_clades()
                if c.name and c.name.split("_")[0] in eskapee_genera]
print(f"ESKAPEE taxa in tree: {len(eskapee_tips)}")

lines = ["DATASET_COLORSTRIP","SEPARATOR TAB",
         "DATASET_LABEL\tESKAPEE","COLOR\t#FF0000",
         "LEGEND_TITLE\tESKAPEE pathogens",
         "LEGEND_SHAPES\t1","LEGEND_COLORS\t#FF0000",
         "LEGEND_LABELS\tESKAPEE pathogen","DATA"]
for tip in eskapee_tips:
    lines.append(f"{tip}\t#FF0000\tESKAPEE")
with open(PATHS["itol_eskapee"],"w") as f: f.write("\n".join(lines))
print(f"✓ {PATHS['itol_eskapee']}")

print("\n══ All files saved to Drive ══")
print("Upload to iTOL:")
print("  1. taninim_silva_tree.treefile — tree file")
print("  2. iTOL_phylum_colors.txt      — phylum ring")
print("  3. iTOL_substrate.txt          — plastic/water ring")
print("  4. iTOL_ESKAPEE_highlight.txt  — ESKAPEE ring")


Source tracking after Matan’s comment

# ============================================================
# FEAST — Taninim River Plastisphere
# Updated: R2→R1, R3→R2, rarefaction depth = 91,000
# Output: C:/Users/metig/Downloads/FEAST_output/
# ============================================================

# --- Step 1: Install packages ------------------------------------------------
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("cozygene/FEAST", dependencies = TRUE, quiet = TRUE)

for (pkg in c("ggplot2", "dplyr", "tidyr", "readr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(FEAST)
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
cat("Libraries loaded\n")

# --- Step 2: Load data -------------------------------------------------------
counts_raw <- read_csv(
  "C:/Users/metig/Downloads/Lig_wf-metagenomics-counts-species.csv",
  show_col_types = FALSE)

meta <- read.delim(
  "C:/Users/metig/Downloads/Ligation_all_metadata (1) (1).txt",
  header = TRUE, stringsAsFactors = FALSE)

# --- Step 3: Fix column names ------------------------------------------------
colnames(meta) <- gsub("^#", "", colnames(meta))
colnames(meta) <- gsub("^X\\.", "", colnames(meta))
meta$SampleID  <- trimws(meta$SampleID)
cat(sprintf("Metadata columns: %s\n", paste(colnames(meta), collapse = ", ")))

# --- Step 4: Build count matrix ----------------------------------------------
tax_cols <- intersect(
  c("total","kingdom","phylum","class","order",
    "family","genus","tax","superkingdom"),
  colnames(counts_raw)
)

count_mat_raw <- counts_raw %>%
  select(-all_of(tax_cols)) %>%
  column_to_rownames(var = "species") %>%
  as.matrix()

# Fix: S1P1 -> S1P4
if ("S1P1" %in% colnames(count_mat_raw)) {
  colnames(count_mat_raw)[colnames(count_mat_raw) == "S1P1"] <- "S1P4"
  cat("Fixed: S1P1 -> S1P4\n")
}

# Transpose: samples x species
count_mat <- t(count_mat_raw)
cat(sprintf("count_mat dimensions: %d samples x %d species\n",
            nrow(count_mat), ncol(count_mat)))
cat("Row names:\n")
print(rownames(count_mat))

# --- Step 5: Filter species --------------------------------------------------
sites <- list(
  R2 = c("R2P1","R2P2","R2P3","R2W1","R2W2","R2W3"),
  R3 = c("R3P1","R3P2","R3P3","R3W1","R3W2","R3W3"),
  S1 = c("S1P4","S1P2","S1P3","S1W1","S1W2","S1W3"),
  Sc = c("ScP1","ScP2","ScP3","ScW1","ScW2","ScW3")
)

pass_total <- colSums(count_mat) >= 50

pass_site <- apply(count_mat, 2, function(sp_counts) {
  any(sapply(sites, function(s) {
    valid <- s[s %in% rownames(count_mat)]
    sum(sp_counts[valid] > 0) >= 2
  }))
})

count_mat <- count_mat[, pass_total & pass_site]
cat(sprintf("After filtering: %d species retained\n", ncol(count_mat)))

# --- Step 6: Verify reads per sample before rarefaction ---------------------
sample_totals <- rowSums(count_mat)
cat("\nReads per sample:\n")
print(sort(sample_totals))
cat(sprintf("\nMinimum reads: %d\n", min(sample_totals)))
cat("Rarefaction depth: 91,000\n")

if (min(sample_totals) < 91000) {
  warning(sprintf(
    "Sample '%s' has only %d reads — below rarefaction depth of 91,000.",
    names(which.min(sample_totals)),
    min(sample_totals)
  ))
}

# --- Step 7: Subsampling function — depth = 91,000 --------------------------
subsample_counts <- function(mat, depth = 91000) {
  t(apply(mat, 1, function(x) {
    if (sum(x) <= depth) return(x)
    as.integer(rmultinom(1, depth, x))
  }))
}

# --- Step 8: FEAST runner function -------------------------------------------
run_feast <- function(meta_feast, label, outdir) {
  samps <- rownames(meta_feast)

  missing <- samps[!samps %in% rownames(count_mat)]
  if (length(missing) > 0) {
    cat(sprintf("WARNING: missing samples: %s\n",
                paste(missing, collapse = ", ")))
    return(NULL)
  }

  mat <- count_mat[samps, ]
  mat[is.na(mat)] <- 0
  mat <- mat[, colSums(mat) > 0]
  mat <- round(mat)
  storage.mode(mat) <- "integer"

  set.seed(42)
  mat <- subsample_counts(mat, depth = 91000)
  colnames(mat) <- colnames(
    count_mat[samps, colSums(count_mat[samps, ]) > 0]
  )

  # Normalize path and create directory
  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!dir.exists(outdir)) {
    stop(sprintf("Could not create output directory: %s", outdir))
  }

  cat(sprintf("Running %s | %d samples x %d species\n",
              label, nrow(mat), ncol(mat)))
  cat(sprintf("Output directory: %s\n", outdir))

  FEAST(C                      = mat,
        metadata               = meta_feast,
        different_sources_flag = 0,
        dir_path               = outdir,
        outfile                = label)
}

# --- Step 9: Create main output folder ---------------------------------------
dir.create("C:/Users/metig/Downloads/FEAST_output",
           showWarnings = FALSE, recursive = TRUE)
cat("Output folder exists:",
    dir.exists("C:/Users/metig/Downloads/FEAST_output"), "\n")

# --- Step 10: Define FEAST metadata ------------------------------------------
# NOTE: Sample IDs in row.names stay as original file names (R2P1 etc.)
# Display labels already use new names (R1, R2)

# Scenario 1: R2 Estuary Plastic (sink)
# Sources: R1 Plastisphere + R1 Water + R2 Water
feast_meta_sc1 <- data.frame(
  Env        = c(rep("R1 Plastisphere", 3),
                 rep("R1 Water", 3),
                 rep("R2 Water", 3),
                 rep("R2 Plastisphere", 3)),
  SourceSink = c(rep("Source", 9), rep("Sink", 3)),
  id         = c(1,1,1, 1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R2W1","R2W2","R2W3",
                 "R3W1","R3W2","R3W3",
                 "R3P1","R3P2","R3P3"),
  stringsAsFactors = FALSE
)

# Scenario 2: R2 Estuary Water (sink)
# Sources: R1 Plastisphere + R1 Water
feast_meta_sc2 <- data.frame(
  Env        = c(rep("R1 Plastisphere", 3),
                 rep("R1 Water", 3),
                 rep("R2 Water", 3)),
  SourceSink = c(rep("Source", 6), rep("Sink", 3)),
  id         = c(1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R2W1","R2W2","R2W3",
                 "R3W1","R3W2","R3W3"),
  stringsAsFactors = FALSE
)

# Scenario 3: S1 Sea Plastic (sink)
# Sources: R1P + R1W + R2P + R2W + S1W
feast_meta_sc3 <- data.frame(
  Env        = c(rep("R1 Plastisphere", 3),
                 rep("R1 Water", 3),
                 rep("R2 Plastisphere", 3),
                 rep("R2 Water", 3),
                 rep("S1 Water", 3),
                 rep("S1 Plastisphere", 3)),
  SourceSink = c(rep("Source", 15), rep("Sink", 3)),
  id         = c(1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R2W1","R2W2","R2W3",
                 "R3P1","R3P2","R3P3",
                 "R3W1","R3W2","R3W3",
                 "S1W1","S1W2","S1W3",
                 "S1P4","S1P2","S1P3"),
  stringsAsFactors = FALSE
)

# Scenario 4: S1 Sea Water (sink)
# Sources: R1P + R2P + R2W
feast_meta_sc4 <- data.frame(
  Env        = c(rep("R1 Plastisphere", 3),
                 rep("R2 Plastisphere", 3),
                 rep("R2 Water", 3),
                 rep("S1 Water", 3)),
  SourceSink = c(rep("Source", 9), rep("Sink", 3)),
  id         = c(1,1,1, 1,1,1, 1,1,1, 1,2,3),
  row.names  = c("R2P1","R2P2","R2P3",
                 "R3P1","R3P2","R3P3",
                 "R3W1","R3W2","R3W3",
                 "S1W1","S1W2","S1W3"),
  stringsAsFactors = FALSE
)

# --- Step 11: Run all 4 scenarios --------------------------------------------
sc1 <- run_feast(feast_meta_sc1, "sc1",
                 "C:/Users/metig/Downloads/FEAST_output/sc1/")
sc2 <- run_feast(feast_meta_sc2, "sc2",
                 "C:/Users/metig/Downloads/FEAST_output/sc2/")
sc3 <- run_feast(feast_meta_sc3, "sc3",
                 "C:/Users/metig/Downloads/FEAST_output/sc3/")
sc4 <- run_feast(feast_meta_sc4, "sc4",
                 "C:/Users/metig/Downloads/FEAST_output/sc4/")

cat("\n=== Scenario 1: R2 Plastisphere (Estuary) ===\n")
print(round(sc1$data_prop * 100, 2))

cat("\n=== Scenario 2: R2 Water (Estuary) ===\n")
print(round(sc2$data_prop * 100, 2))

cat("\n=== Scenario 3: S1 Plastisphere (Sea) ===\n")
print(round(sc3$data_prop * 100, 2))

cat("\n=== Scenario 4: S1 Water (Sea) ===\n")
print(round(sc4$data_prop * 100, 2))

# --- Step 12: Summarize results ----------------------------------------------
summarize_scenario <- function(data_prop, sink_label) {
  df <- as.data.frame(data_prop[, "pred_emnoise_all", drop = FALSE])
  df$source_replicate <- rownames(df)
  df$Env <- gsub("_[^_]+$", "", df$source_replicate)
  env_sums <- df %>%
    group_by(Env) %>%
    summarise(Proportion = sum(pred_emnoise_all) * 100, .groups = "drop")
  env_sums$Proportion <- env_sums$Proportion /
    sum(env_sums$Proportion) * 100
  env_sums$Sink <- sink_label
  return(env_sums)
}

sc1_sum <- summarize_scenario(sc1$data_prop, "R2 Plastisphere (Estuary)")
sc2_sum <- summarize_scenario(sc2$data_prop, "R2 Water (Estuary)")
sc3_sum <- summarize_scenario(sc3$data_prop, "S1 Plastisphere (Sea)")
sc4_sum <- summarize_scenario(sc4$data_prop, "S1 Water (Sea)")

plot_df <- bind_rows(sc1_sum, sc2_sum, sc3_sum, sc4_sum)

# --- Step 13: Plot -----------------------------------------------------------
plot_df$Sink <- factor(plot_df$Sink, levels = c(
  "R2 Plastisphere (Estuary)",
  "R2 Water (Estuary)",
  "S1 Plastisphere (Sea)",
  "S1 Water (Sea)"
))

source_order <- c(
  "R1 Plastisphere", "R2 Plastisphere",
  "R1 Water",        "R2 Water",
  "S1 Water",        "unknown"
)
plot_df$Env <- factor(plot_df$Env, levels = source_order)

plot_df$Label <- ifelse(
  plot_df$Proportion >= 3,
  paste0(round(plot_df$Proportion, 1), "%"), ""
)

source_colors <- c(
  "R1 Plastisphere" = "#7B2D8B",
  "R2 Plastisphere" = "#2B5EA7",
  "R1 Water"        = "#E8638C",
  "R2 Water"        = "#6B1A1A",
  "S1 Water"        = "#148F77",
  "unknown"         = "#D4860A"
)

p_st <- ggplot(plot_df,
               aes(x = Sink, y = Proportion, fill = Env)) +
  geom_bar(stat = "identity", width = 0.55,
           color = "black", linewidth = 0.8) +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = source_colors, name = "Source") +
  scale_x_discrete(limits = rev(c(
    "R2 Plastisphere (Estuary)",
    "R2 Water (Estuary)",
    "S1 Plastisphere (Sea)",
    "S1 Water (Sea)"
  ))) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 101),
    breaks = seq(0, 100, 20),
    labels = paste0(seq(0, 100, 20), "%")
  ) +
  coord_flip() +
  labs(
    title    = "Microbial Source Tracking — Taninim River to Mediterranean Sea",
    subtitle = paste0(
      "FEAST Bayesian source tracking | Rarefied to 91,000 reads\n",
      "R1: Upstream Taninim River | R2: Taninim River Estuary | ",
      "S1: Connected EMS Beach | Sc: Caesarea Recreational Beach"
    ),
    x = NULL,
    y = "Source Contribution (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 9, color = "grey40"),
    axis.text.y       = element_text(size = 11, color = "black",
                                     face = "bold"),
    axis.text.x       = element_text(size = 11, color = "black",
                                     face = "bold"),
    axis.title.x      = element_text(size = 12, color = "black",
                                     face = "bold"),
    axis.line         = element_line(color = "black", linewidth = 1.2),
    axis.ticks        = element_line(color = "black", linewidth = 1.0),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title      = element_text(face = "bold", size = 11),
    legend.text       = element_text(size = 10),
    legend.position   = "right",
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

print(p_st)

# --- Step 14: Save outputs ---------------------------------------------------
ggsave(
  "C:/Users/metig/Downloads/FEAST_output/FEAST_sourcetracking_4scenarios.png",
  p_st, width = 10, height = 6, dpi = 300, bg = "white"
)

write.csv(
  plot_df,
  "C:/Users/metig/Downloads/FEAST_output/FEAST_summary.csv",
  row.names = FALSE
)

cat("Saved: C:/Users/metig/Downloads/FEAST_output/FEAST_sourcetracking_4scenarios.png\n")
cat("Saved: C:/Users/metig/Downloads/FEAST_output/FEAST_summary.csv\n")
cat("\nDone!\n")


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

tree <- read.tree("D:/taninim_tree2.treefile")
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
           "D:/itol_ring_substrate.txt")
cat("✅ D:/itol_ring_substrate.txt\n")


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
           "D:/itol_ring_site.txt")
cat("✅ D:/itol_ring_site.txt\n")


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

writeLines(branch_lines, "D:/itol_branches_thick.txt")
cat("✅ D:/itol_branches_thick.txt\n")


# ── STEP 5: Dominant family ring ─────────────────────────────
tax <- read.csv("D:/Ligation_abundance_table_counts-species.csv",
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

writeLines(ring_lines, "D:/itol_family_ring.txt")
cat("✅ D:/itol_family_ring.txt\n")


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
writeLines(c(legend_lines, data_lines), "D:/itol_site_legend.txt")
cat("✅ D:/itol_site_legend.txt\n")


# ── SUMMARY ──────────────────────────────────────────────────
cat("\n============================================================\n")
cat("Files generated — drag onto iTOL in this order:\n")
cat("  1. D:/itol_branches_thick.txt   → branch colors by site\n")
cat("  2. D:/itol_ring_substrate.txt   → inner ring: Plastic/Water\n")
cat("  3. D:/itol_family_ring.txt      → outer ring: dominant family\n")
cat("  4. D:/itol_site_legend.txt      → legend for branch colors\n")
cat("============================================================\n")


Taninim Plastisphere Phylogenetic Tree - Run R Scripts
# Run this in PowerShell to generate all iTOL annotation files
# ============================================================

# Save the R script first, then run it with Rscript
$rScript = @"
library(ape)
library(phangorn)

tree <- read.tree("D:/taninim_tree2.treefile")
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
           "D:/itol_ring_substrate.txt")
cat("OK D:/itol_ring_substrate.txt\n")

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
           "D:/itol_ring_site.txt")
cat("OK D:/itol_ring_site.txt\n")

# ── STEP 4: Branch colors ────────────────────────────────────
branch_lines <- c("TREE_COLORS","SEPARATOR TAB","","DATA")
for (i in seq_along(tips)) {
  branch_lines <- c(branch_lines,
    paste(tips[i], "branch", site_colors[df\$site[i]], "normal", "4", sep="\t"))
}
writeLines(branch_lines, "D:/itol_branches_thick.txt")
cat("OK D:/itol_branches_thick.txt\n")

# ── STEP 5: Dominant family ring ─────────────────────────────
tax <- read.csv("D:/Ligation_abundance_table_counts-species.csv",
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
writeLines(ring_lines, "D:/itol_family_ring.txt")
cat("OK D:/itol_family_ring.txt\n")

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
           "D:/itol_site_legend.txt")
cat("OK D:/itol_site_legend.txt\n")

cat("\nAll files generated successfully!\n")
cat("Drag onto iTOL in this order:\n")
cat("  1. itol_branches_thick.txt\n")
cat("  2. itol_ring_substrate.txt\n")
cat("  3. itol_family_ring.txt\n")
cat("  4. itol_site_legend.txt\n")
"@

# Save R script to temp file
$rScript | Out-File -FilePath "D:\itol_generate.R" -Encoding UTF8

# Run with Rscript
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" "D:\itol_generate.R"

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

