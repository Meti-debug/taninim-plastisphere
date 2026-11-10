# Non-marine taxa in river-mouth water (Results 3.3, Supplementary-6)
#
# Taxa detected (>= 1 read) in at least two of three S1 water replicates and in
# at least two of three R2 plastic or R2 water replicates, and with no reads in
# any Caesarea water replicate. Reports their number, their summed share of S1
# water reads, and the twelve most abundant.
#
# Uses the FEAST-filtered table rarefied to 91,000 reads (seed 42), as in
# r/09_feast_sensitivity.R. Sequencing labels: R2* = R1, R3* = R2.
#
# Run from the repository root:
#   Rscript r/10_river_signature.R

source("r/00_paths.R")
DEPTH <- 91000

counts_raw <- read.csv(file.path(DATA, "counts.csv"), check.names = FALSE)
tax_cols <- intersect(c("total", "kingdom", "phylum", "class", "order", "family",
                        "genus", "tax", "superkingdom"), colnames(counts_raw))
mat <- as.matrix(counts_raw[, setdiff(colnames(counts_raw), c(tax_cols, "species"))])
rownames(mat) <- counts_raw$species
if ("S1P1" %in% colnames(mat)) colnames(mat)[colnames(mat) == "S1P1"] <- "S1P4"
mat[is.na(mat)] <- 0
X <- t(mat)

sites <- list(R2 = c("R2P1","R2P2","R2P3","R2W1","R2W2","R2W3"),
              R3 = c("R3P1","R3P2","R3P3","R3W1","R3W2","R3W3"),
              S1 = c("S1P4","S1P2","S1P3","S1W1","S1W2","S1W3"),
              Sc = c("ScP1","ScP2","ScP3","ScW1","ScW2","ScW3"))
keep <- colSums(X) >= 50 & apply(X, 2, function(x) any(sapply(sites, function(s) sum(x[s] > 0) >= 2)))
X <- X[, keep]
set.seed(42)
storage.mode(X) <- "integer"
R <- t(apply(X, 1, function(x) if (sum(x) <= DEPTH) x else as.integer(rmultinom(1, DEPTH, x))))
colnames(R) <- colnames(X)
cat("Species after filtering:", ncol(R), "\n")

lib <- function(p) paste0(p, 1:3)
P <- function(g) R[lib(g), , drop = FALSE] / DEPTH
n_det <- function(g) colSums(R[lib(g), , drop = FALSE] > 0)

# R3* = manuscript R2
estuary <- n_det("R3P") >= 2 | n_det("R3W") >= 2
in_s1w  <- n_det("S1W") >= 2
no_scw  <- colSums(R[lib("ScW"), , drop = FALSE]) == 0
set <- estuary & in_s1w & no_scw

s1w <- colMeans(P("S1W"))
cat(sprintf("Taxa: %d; share of S1 water reads: %.2f%%\n", sum(set), 100 * sum(s1w[set])))
cat("Per S1 water replicate (%):",
    sprintf("%.2f", 100 * rowSums(P("S1W")[, set, drop = FALSE])), "\n\n")

top <- names(sort(s1w[set], decreasing = TRUE))[1:12]
print(data.frame(species    = top,
                 S1W_pct    = round(100 * s1w[top], 3),
                 R2P_pct    = round(100 * colMeans(P("R3P"))[top], 3),
                 R2W_pct    = round(100 * colMeans(P("R3W"))[top], 3)),
      row.names = FALSE)
