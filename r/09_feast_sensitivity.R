# FEAST sensitivity analysis for the river-mouth water sink
#
# Re-runs the S1 water source-tracking scenario with the source pool changed:
#   A  baseline      R1 plastic, R2 plastic, R2 water (as reported)
#   B  drop R2P      R1 plastic, R2 water
#   C  drop R1P      R2 plastic, R2 water
#   D  + marine      baseline plus Caesarea water as a marine endmember
#   E  full river    R1 plastic, R1 water, R2 plastic, R2 water
#   F  no plastic    R2 water
#   G  river water   R1 water, R2 water
#   H  water+marine  R2 water, Caesarea water
#
# If R2 plastic carries a signal that no other source shares, removing it (B)
# should move most of its share into Unknown. If R1 plastic or R2 water absorb
# it instead, the attribution rests on taxa shared between sources.
#
# Species filtering, rarefaction to 91,000 reads (seed 42) and the multiple-sink
# configuration with shared sources follow Methods 2.14. Each scenario is run
# five times with seeds 1 to 5. Contributions are summed over the source
# replicates of an environment, then summarised over 3 sinks x 5 runs.
#
# Sample names in counts.csv use the sequencing labels:
#   R2* = manuscript R1, R3* = manuscript R2, S1P4 = manuscript S1P1.
#
# Run from the repository root:
#   Rscript r/09_feast_sensitivity.R

source("r/00_paths.R")

library(FEAST)

counts_file <- file.path(DATA, "counts.csv")
out_dir     <- normalizePath(file.path(OUT, "feast_sensitivity"), mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DEPTH <- 91000
RUNS  <- 1:5

# -- 1. Count matrix, samples x species ----------------------------------------
counts_raw <- read.csv(counts_file, check.names = FALSE)
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
cat(sprintf("Species retained: %d; reads retained: %.2f%%\n", ncol(count_mat),
            100 * sum(count_mat) / sum(mat)))

# -- 3. Rarefy every library to 91,000 reads -----------------------------------
set.seed(42)
storage.mode(count_mat) <- "integer"
rar <- t(apply(count_mat, 1, function(x)
  if (sum(x) <= DEPTH) x else as.integer(rmultinom(1, DEPTH, x))))
colnames(rar) <- colnames(count_mat)
storage.mode(rar) <- "integer"
cat("Library depths after rarefaction:", range(rowSums(rar)), "\n")

# -- 4. Scenarios --------------------------------------------------------------
lib <- function(prefix) paste0(prefix, 1:3)
src <- list(
  "R1 plastic" = lib("R2P"), "R1 water" = lib("R2W"),
  "R2 plastic" = lib("R3P"), "R2 water" = lib("R3W"),
  "Sc water"   = lib("ScW"))
sink <- lib("S1W")

scenarios <- list(
  A_baseline   = c("R1 plastic", "R2 plastic", "R2 water"),
  B_drop_R2P   = c("R1 plastic", "R2 water"),
  C_drop_R1P   = c("R2 plastic", "R2 water"),
  D_add_marine = c("R1 plastic", "R2 plastic", "R2 water", "Sc water"),
  E_full_river = c("R1 plastic", "R1 water", "R2 plastic", "R2 water"),
  F_water_only = c("R2 water"),
  G_river_water = c("R1 water", "R2 water"),
  H_water_marine = c("R2 water", "Sc water"))

run_one <- function(envs, seed, label, scen) {
  ids  <- c(unlist(src[envs]), sink)
  meta <- data.frame(
    Env        = c(rep(envs, each = 3), rep("S1 water", 3)),
    SourceSink = c(rep("Source", 3 * length(envs)), rep("Sink", 3)),
    id         = c(rep(NA, 3 * length(envs)), 1:3),
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
  agg <- sapply(c(envs, "Unknown"), function(e)
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
write.csv(summ, file.path(out_dir, "feast_sensitivity_summary.csv"), row.names = FALSE)
cat("\nWritten:", file.path(out_dir, "feast_sensitivity_summary.csv"), "\n")
