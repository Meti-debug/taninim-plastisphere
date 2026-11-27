# FEAST source tracking: the four primary source-sink scenarios (Figure 4E)
#
#   P1  R2 plastic sink   sources R1 plastic, R1 water, R2 water
#   P2  R2 water sink     sources R1 plastic, R1 water
#   P3  S1 plastic sink   sources R1 plastic, R1 water, R2 plastic, R2 water, S1 water
#   P4  S1 water sink     sources R1 plastic, R2 plastic, R2 water
#
# Species filtering, rarefaction to 91,000 reads (seed 42) and the multiple-sink
# configuration with shared sources are the same as in sensitivity_analysis.R,
# so P4 is identical to its baseline scenario A. Each scenario is run five
# times with seeds 1 to 5. Contributions are summed over the source replicates
# of an environment, then summarised as mean and SD over 3 sinks x 5 runs.
#
# Sample names are handled as sequencing labels (see analysis/paths.R):
#   R2* = manuscript R1, R3* = manuscript R2, S1P4 = manuscript S1P1.
#
# Run from the repository root:
#   Rscript analysis/06_FEAST/primary_analysis.R

source("analysis/paths.R")
OUT <- step_output("06_FEAST")

library(FEAST)
library(ggplot2)

counts_file <- file.path(DATA, "counts.csv")
out_dir     <- normalizePath(file.path(OUT, "feast_primary"), mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DEPTH <- 91000
RUNS  <- 1:5

# -- 1. Count matrix, samples x species ----------------------------------------
counts_raw <- seq_labels(read.csv(counts_file, check.names = FALSE))
tax_cols <- intersect(
  c("total", "kingdom", "phylum", "class", "order", "family", "genus", "tax",
    "superkingdom"), colnames(counts_raw))
mat <- as.matrix(counts_raw[, setdiff(colnames(counts_raw), c(tax_cols, "species"))])
rownames(mat) <- counts_raw$species
if ("S1P1" %in% colnames(mat)) colnames(mat)[colnames(mat) == "S1P1"] <- "S1P4"
mat[is.na(mat)] <- 0
count_mat <- t(mat)

# -- 2. Species filter: >= 50 reads in total and >= 2 samples at one site ------
sites <- list(
  R2 = c("R2P1", "R2P2", "R2P3", "R2W1", "R2W2", "R2W3"),
  R3 = c("R3P1", "R3P2", "R3P3", "R3W1", "R3W2", "R3W3"),
  S1 = c("S1P4", "S1P2", "S1P3", "S1W1", "S1W2", "S1W3"),
  Sc = c("ScP1", "ScP2", "ScP3", "ScW1", "ScW2", "ScW3"))
pass_total <- colSums(count_mat) >= 50
pass_site  <- apply(count_mat, 2, function(x)
  any(sapply(sites, function(s) sum(x[s] > 0) >= 2)))
count_mat <- count_mat[, pass_total & pass_site]
cat(sprintf("Species retained: %d of %d; reads retained: %.2f%% of %s\n",
            ncol(count_mat), ncol(t(mat)), 100 * sum(count_mat) / sum(mat),
            format(sum(mat), big.mark = ",")))

# -- 3. Rarefy every library to 91,000 reads -----------------------------------
set.seed(42)
storage.mode(count_mat) <- "integer"
rar <- t(apply(count_mat, 1, function(x)
  if (sum(x) <= DEPTH) x else as.integer(rmultinom(1, DEPTH, x))))
colnames(rar) <- colnames(count_mat)
storage.mode(rar) <- "integer"
cat("Library depths after rarefaction:", range(rowSums(rar)), "\n")

# -- 4. Scenarios --------------------------------------------------------------
lib <- function(prefix) if (prefix == "S1P") c("S1P4", "S1P2", "S1P3") else paste0(prefix, 1:3)
envlib <- c("R1 plastic" = "R2P", "R1 water" = "R2W", "R2 plastic" = "R3P",
            "R2 water" = "R3W", "S1 plastic" = "S1P", "S1 water" = "S1W")

scenarios <- list(
  P1_R2P_sink = list(sink = "R2 plastic", src = c("R1 plastic", "R1 water", "R2 water")),
  P2_R2W_sink = list(sink = "R2 water",   src = c("R1 plastic", "R1 water")),
  P3_S1P_sink = list(sink = "S1 plastic", src = c("R1 plastic", "R1 water", "R2 plastic",
                                                  "R2 water", "S1 water")),
  P4_S1W_sink = list(sink = "S1 water",   src = c("R1 plastic", "R2 plastic", "R2 water")))

run_one <- function(sc, seed, label, scen) {
  ids  <- c(unlist(lapply(envlib[sc$src], lib)), lib(envlib[[sc$sink]]))
  meta <- data.frame(
    Env        = c(rep(sc$src, each = 3), rep(sc$sink, 3)),
    SourceSink = c(rep("Source", 3 * length(sc$src)), rep("Sink", 3)),
    id         = c(rep(NA, 3 * length(sc$src)), 1:3),
    row.names  = ids, stringsAsFactors = FALSE)
  set.seed(seed)
  wd <- getwd()                         # FEAST calls setwd(dir_path)
  FEAST(C = rar[ids, ], metadata = meta, different_sources_flag = 0,
        dir_path = out_dir, outfile = label)
  setwd(wd)
  P <- read.delim(file.path(out_dir, paste0(label, "_source_contributions_matrix.txt")),
                  check.names = FALSE)
  # columns are "<sample>_<Env>" plus "Unknown"; sum the source replicates
  env_of <- sub("^[^_]+_", "", colnames(P))
  agg <- sapply(c(sc$src, "Unknown"), function(e)
    rowSums(P[, env_of == e, drop = FALSE], na.rm = TRUE))
  data.frame(scenario = scen, seed = seed, sink = rownames(P), agg * 100,
             check.names = FALSE, row.names = NULL)
}

res <- list()
for (sc in names(scenarios)) for (s in RUNS) {
  cat(sprintf("\n== %s, seed %d ==\n", sc, s))
  res[[paste(sc, s)]] <- run_one(scenarios[[sc]], s, paste0(sc, "_seed", s), sc)
}

# -- 5. Summary: mean and SD over 3 sinks x 5 runs -----------------------------
long <- do.call(rbind, lapply(res, function(r) {
  v <- r[, -(1:3), drop = FALSE]
  data.frame(scenario = r$scenario[1], source = rep(colnames(v), each = nrow(v)),
             value = unlist(v), row.names = NULL)
}))
summ <- aggregate(value ~ scenario + source, data = long,
                  FUN = function(x) c(mean = mean(x), sd = sd(x)))
summ <- data.frame(scenario = summ$scenario, source = summ$source,
                   mean = round(summ$value[, "mean"], 1),
                   sd = round(summ$value[, "sd"], 1))
summ <- summ[order(summ$scenario, -summ$mean), ]
print(summ, row.names = FALSE)
write.csv(summ, file.path(out_dir, "feast_primary_summary.csv"), row.names = FALSE)
write.csv(long, file.path(out_dir, "feast_primary_all_runs.csv"), row.names = FALSE)

# -- 6. Stacked bars, one per sink (the values drawn in Figure 4E) -------------
sink_lab <- c(P1_R2P_sink = "R2P", P2_R2W_sink = "R2W",
              P3_S1P_sink = "S1P", P4_S1W_sink = "S1W")
src_lab  <- c("R1 plastic" = "R1P", "R1 water" = "R1W", "R2 plastic" = "R2P",
              "R2 water" = "R2W", "S1 water" = "S1W", "Unknown" = "Unknown")
pd <- data.frame(sink   = factor(sink_lab[summ$scenario], levels = rev(sink_lab)),
                 source = factor(src_lab[summ$source], levels = src_lab),
                 mean   = summ$mean, sd = summ$sd)
pd$label <- ifelse(pd$mean >= 10, sprintf("%.1f \u00b1 %.1f%%", pd$mean, pd$sd), "")
p <- ggplot(pd, aes(x = sink, y = mean, fill = source)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label), size = 2.6,
            position = position_stack(vjust = 0.5, reverse = TRUE)) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = c(R1P = "#F28E2B", R1W = "#E5C3E8", R2P = "#4E79A7",
                               R2W = "#A9C0DA", S1W = "#A8D5A2", Unknown = "#BDBDBD"),
                    name = NULL) +
  labs(x = NULL, y = "Source contribution (%)") +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))
ggsave(file.path(out_dir, "Figure4E_FEAST.pdf"), p, width = 6, height = 3)
ggsave(file.path(out_dir, "Figure4E_FEAST.png"), p, width = 6, height = 3, dpi = 300)

cat("\nWritten:", file.path(out_dir, "feast_primary_summary.csv"), "\n")
cat("        ", file.path(out_dir, "Figure4E_FEAST.pdf/.png"), "\n")
