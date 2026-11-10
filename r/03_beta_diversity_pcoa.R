# Beta diversity pcoa
#
# Bray-Curtis dissimilarity, PCoA ordination and pairwise dissimilarity panels.
#
# Extracted from r/figures_as_run.R, lines 1436-1675, which is the session code
# deposited verbatim. Paths were made relative; nothing else was changed.
#
# Run from the repository root:
#   Rscript r/03_beta_diversity_pcoa.R

source("r/00_paths.R")

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── 1. Load data ──────────────────────────────────────────────────────────────
meta <- read.table(
  file.path(DATA, "metadata.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sp_tab <- read.csv(
  file.path(DATA, "counts.csv"),
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

