# River signature at the river mouth versus the reference beach
#
# Compares Bray-Curtis dissimilarity between river samples and the two marine
# sites: S1 (river mouth) and Sc (Caesarea, not affected by river discharge).
# If plastic carried river-derived communities to the mouth, river plastic would
# be closer to S1 plastic than to Sc plastic; under selection by the surface the
# two would be equally similar. The same comparison is made for water.
# Statistic: mean distance from river samples to Sc minus mean distance to S1,
# with an exact permutation of the six marine labels (20 splits).
# Reported in Results 3.2.
#
# Run from the repository root:
#   Rscript r/13_marine_site_comparison.R

source("r/00_paths.R")
suppressMessages(library(vegan))

x <- read.csv(file.path(DATA, "counts.csv"), check.names = FALSE)
tax <- intersect(c("total", "kingdom", "phylum", "class", "order", "family",
                   "genus", "tax", "superkingdom", "species"), colnames(x))
m <- as.matrix(x[, setdiff(colnames(x), tax)]); m[is.na(m)] <- 0
D <- as.matrix(vegdist(t(m), "bray"))           # read counts, as in the PERMANOVA

# sequencing labels: R2* = R1, R3* = R2, S1P4 = S1P1
lib <- function(p) if (p == "S1P") c("S1P4", "S1P2", "S1P3") else paste0(p, 1:3)
md  <- function(a, b) median(D[lib(a), lib(b)])

compare <- function(river, s1, sc, label) {
  for (r in river) cat(sprintf("  %s: to S1 %.3f | to Sc %.3f\n", r, md(r, s1), md(r, sc)))
  riv <- unlist(lapply(river, lib)); mar <- c(lib(s1), lib(sc))
  stat <- function(s1_ids) mean(D[riv, setdiff(mar, s1_ids)]) - mean(D[riv, s1_ids])
  obs  <- stat(lib(s1))
  null <- combn(mar, 3, stat)
  cat(sprintf("  %s closer to S1 than to Sc by %.3f; exact permutation p = %.2f (%d splits)\n\n",
              label, obs, mean(null >= obs), length(null)))
}
cat("Plastic (R2P = R1 plastic, R3P = R2 plastic)\n")
compare(c("R2P", "R3P"), "S1P", "ScP", "river plastic")
cat("Water (R2W = R1 water, R3W = R2 water)\n")
compare(c("R2W", "R3W"), "S1W", "ScW", "river water")
cat(sprintf("S1 plastic vs Sc plastic: %.3f\n", md("S1P", "ScP")))
