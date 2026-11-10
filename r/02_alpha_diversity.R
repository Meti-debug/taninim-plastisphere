# Alpha diversity
#
# Shannon diversity per sample type, with Wilcoxon tests and boxplots.
#
# Extracted from r/figures_as_run.R, lines 1676-1860, which is the session code
# deposited verbatim. Paths were made relative; nothing else was changed.
#
# Run from the repository root:
#   Rscript r/02_alpha_diversity.R

source("r/00_paths.R")

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

