# FEAST source tracking at genus level
#
# Repeats the four primary source-sink scenarios and the eight S1 water source
# pools of r/09_feast_sensitivity.R on genus-level counts, to test whether the
# attributions depend on species-level assignment. Species counts are summed by
# genus; filtering, rarefaction to 91,000 reads (seed 42), the multiple-sink
# configuration and seeds 1 to 5 are as at species level.
#
# Sequencing labels: R2* = manuscript R1, R3* = manuscript R2, S1P4 = S1P1.
#
# Run from the repository root:
#   Rscript r/11_feast_genus.R

source("r/00_paths.R")
library(FEAST)

out_dir <- normalizePath(file.path(OUT, "feast_genus"), mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
DEPTH <- 91000
RUNS  <- 1:5

# -- 1. Genus x sample counts ---------------------------------------------------
counts_raw <- read.csv(file.path(DATA, "counts.csv"), check.names = FALSE)
tax_cols <- intersect(c("total", "kingdom", "phylum", "class", "order", "family",
                        "genus", "tax", "superkingdom", "species"), colnames(counts_raw))
mat <- as.matrix(counts_raw[, setdiff(colnames(counts_raw), tax_cols)])
mat[is.na(mat)] <- 0
mat <- rowsum(mat, group = counts_raw$genus)          # sum species within genus
if ("S1P1" %in% colnames(mat)) colnames(mat)[colnames(mat) == "S1P1"] <- "S1P4"
count_mat <- t(mat)

# -- 2. Filter: >= 50 reads in total and >= 2 samples at one site --------------
sites <- list(
  R2 = c("R2P1", "R2P2", "R2P3", "R2W1", "R2W2", "R2W3"),
  R3 = c("R3P1", "R3P2", "R3P3", "R3W1", "R3W2", "R3W3"),
  S1 = c("S1P4", "S1P2", "S1P3", "S1W1", "S1W2", "S1W3"),
  Sc = c("ScP1", "ScP2", "ScP3", "ScW1", "ScW2", "ScW3"))
pass_total <- colSums(count_mat) >= 50
pass_site  <- apply(count_mat, 2, function(x)
  any(sapply(sites, function(s) sum(x[s] > 0) >= 2)))
count_mat <- count_mat[, pass_total & pass_site]
cat(sprintf("Genera retained: %d of %d; reads retained: %.2f%%\n",
            ncol(count_mat), nrow(mat), 100 * sum(count_mat) / sum(mat)))

# -- 3. Rarefy to 91,000 reads -------------------------------------------------
set.seed(42)
storage.mode(count_mat) <- "integer"
rar <- t(apply(count_mat, 1, function(x)
  if (sum(x) <= DEPTH) x else as.integer(rmultinom(1, DEPTH, x))))
colnames(rar) <- colnames(count_mat)
storage.mode(rar) <- "integer"

# -- 4. Scenarios: sink and its sources ----------------------------------------
lib <- function(prefix) if (prefix == "S1P") c("S1P4", "S1P2", "S1P3") else paste0(prefix, 1:3)
envlib <- c("R1 plastic" = "R2P", "R1 water" = "R2W", "R2 plastic" = "R3P",
            "R2 water" = "R3W", "S1 plastic" = "S1P", "S1 water" = "S1W", "Sc water" = "ScW")
scenarios <- list(
  P1_R2P_sink    = list(sink = "R2 plastic", src = c("R1 plastic", "R1 water", "R2 water")),
  P2_R2W_sink    = list(sink = "R2 water",   src = c("R1 plastic", "R1 water")),
  P3_S1P_sink    = list(sink = "S1 plastic", src = c("R1 plastic", "R1 water", "R2 plastic", "R2 water", "S1 water")),
  A_baseline     = list(sink = "S1 water",   src = c("R1 plastic", "R2 plastic", "R2 water")),
  B_drop_R2P     = list(sink = "S1 water",   src = c("R1 plastic", "R2 water")),
  C_drop_R1P     = list(sink = "S1 water",   src = c("R2 plastic", "R2 water")),
  D_add_marine   = list(sink = "S1 water",   src = c("R1 plastic", "R2 plastic", "R2 water", "Sc water")),
  E_full_river   = list(sink = "S1 water",   src = c("R1 plastic", "R1 water", "R2 plastic", "R2 water")),
  F_water_only   = list(sink = "S1 water",   src = c("R2 water")),
  G_river_water  = list(sink = "S1 water",   src = c("R1 water", "R2 water")),
  H_water_marine = list(sink = "S1 water",   src = c("R2 water", "Sc water")))

run_one <- function(sc, seed, label, name) {
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
  env_of <- sub("^[^_]+_", "", colnames(P))
  agg <- sapply(c(sc$src, "Unknown"), function(e)
    rowSums(P[, env_of == e, drop = FALSE], na.rm = TRUE))
  data.frame(scenario = name, seed = seed, sink = rownames(P), agg * 100,
             check.names = FALSE, row.names = NULL)
}

res <- list()
for (nm in names(scenarios)) for (s in RUNS) {
  cat(sprintf("\n== %s, seed %d ==\n", nm, s))
  res[[paste(nm, s)]] <- run_one(scenarios[[nm]], s, paste0(nm, "_seed", s), nm)
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
                   mean = round(summ$value[, "mean"], 1), sd = round(summ$value[, "sd"], 1))
summ <- summ[order(summ$scenario, -summ$mean), ]
print(summ, row.names = FALSE)
write.csv(summ, file.path(out_dir, "feast_genus_summary.csv"), row.names = FALSE)
cat("\nWritten:", file.path(out_dir, "feast_genus_summary.csv"), "\n")
