# Source tracking feast
#
# FEAST source tracking for the four source-sink scenarios, and its stacked bars.
#
# Taken from the R session code used for the paper; paths were made
# relative and nothing else was changed.
#
# Run from the repository root:
#   Rscript analysis/r/06_source_tracking_feast.R

source("analysis/r/00_paths.R")

library(FEAST)
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
cat("Libraries loaded\n")

# --- Step 2: Load data ---
counts_raw <- read_csv(file.path(DATA, "counts.csv"),
                       show_col_types = FALSE)
meta       <- read.delim(file.path(DATA, "metadata.tsv"),
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
sc1 <- run_feast(feast_meta_sc1, "sc1", file.path(OUT, ""))
sc2 <- run_feast(feast_meta_sc2, "sc2", file.path(OUT, ""))
sc3 <- run_feast(feast_meta_sc3, "sc3", file.path(OUT, ""))
sc4 <- run_feast(feast_meta_sc4, "sc4", file.path(OUT, ""))

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
ggsave(file.path(OUT, "FEAST_sourcetracking_4scenarios.png"), p_st,
       width = 10, height = 6, dpi = 300, bg = "white")

write.csv(plot_df, file.path(OUT, "FEAST_summary.csv"), row.names = FALSE)
cat("Saved: output/FEAST_sourcetracking_4scenarios.png\n")
cat("Saved: output/FEAST_summary.csv\n")
cat("\nDone!\n")
