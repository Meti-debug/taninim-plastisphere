# Source tracking feast
#
# FEAST source tracking for the four source-sink scenarios, and its stacked bars.
#
# Extracted from r/figures_as_run.R, lines 2503-3060, which is the session code
# deposited verbatim. Paths were made relative; nothing else was changed.
#
# Run from the repository root:
#   Rscript r/06_source_tracking_feast.R

source("r/00_paths.R")

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

