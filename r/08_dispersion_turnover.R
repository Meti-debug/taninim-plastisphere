# Dispersion and river-sea turnover
#
# Homogeneity of multivariate dispersion (betadisper / permutest) for the
# Bray-Curtis matrix used in the PERMANOVA, and a within-site permutation test
# of the difference in river-sea turnover between plastic-bound and planktonic
# communities, and a check of the river or sea x sample type PERMANOVA.
# Reported in Results 3.2 and Supplementary-1.
#
# The count matrix, site relabelling and Bray-Curtis matrix are built exactly as
# in r/03_beta_diversity_pcoa.R.
#
# Run from the repository root:
#   Rscript r/08_dispersion_turnover.R

source("r/00_paths.R")

library(vegan)
library(dplyr)

set.seed(42)
NPERM <- 9999

# -- 1. Load data and build the Bray-Curtis matrix -----------------------------
meta <- read.table(
  file.path(DATA, "metadata.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sp_tab <- read.csv(
  file.path(DATA, "counts.csv"),
  header = TRUE, check.names = FALSE)

sp_counts <- sp_tab[, meta$SampleID]
rownames(sp_counts) <- sp_tab$species
sp_counts_clean <- sp_counts[rowSums(sp_counts) > 0, ]

meta_clean <- meta
meta_clean$Site <- dplyr::recode(meta_clean$Site, "R2" = "R1", "R3" = "R2")

bray <- vegdist(t(sp_counts_clean), method = "bray")
stopifnot(all(labels(bray) == meta_clean$SampleID))

site <- factor(meta_clean$Site, levels = c("R1", "R2", "S1", "Sc"))
type <- factor(meta_clean$SampleType)
rs   <- factor(ifelse(site %in% c("R1", "R2"), "river", "sea"))
print(table(site, type))

# Check against the manuscript: site R2 0.307, sample type R2 0.217, model R2 0.524.
# This call also fixes the random-number stream, so the p-values below match the
# reported ones.
print(adonis2(bray ~ site + type, permutations = NPERM))

# -- 2. Dispersion: site x sample type (8 groups of 3) -------------------------
cat("\n## Site x sample type\n")
bd_grp <- betadisper(bray, interaction(site, type))
print(permutest(bd_grp, permutations = NPERM))

# -- 3. Dispersion: plastic vs water, all sites --------------------------------
cat("\n## Sample type, all sites\n")
bd_type <- betadisper(bray, type)
print(permutest(bd_type, permutations = NPERM))
print(tapply(bd_type$distances, type, median))

# -- 4. Dispersion: river or sea x sample type (4 groups of 6) -----------------
cat("\n## River or sea x sample type\n")
g_rs  <- interaction(rs, type)
bd_rs <- betadisper(bray, g_rs)
print(permutest(bd_rs, permutations = NPERM, pairwise = TRUE))
print(tapply(bd_rs$distances, g_rs, median))

# -- 5. Dispersion: plastic vs water within each site --------------------------
# Three replicates per sample type, so permutest enumerates all 720 permutations.
for (s in levels(site)) {
  k <- site == s
  cat("\n## Sample type within", s, "\n")
  print(permutest(betadisper(as.dist(as.matrix(bray)[k, k]), droplevels(type[k])),
                  permutations = 999))
}

# -- 6. River-sea turnover: plastic vs water -----------------------------------
# Statistic: median Bray-Curtis dissimilarity of water river-sea pairs minus that
# of plastic river-sea pairs. Sample-type labels are permuted within site.
set.seed(42)
D <- as.matrix(bray)
cross <- outer(rs == "river", rs == "sea")

turnover_diff <- function(tp) {
  w <- cross & outer(tp == "Water",   tp == "Water")
  p <- cross & outer(tp == "Plastic", tp == "Plastic")
  median(D[w]) - median(D[p])
}

obs <- turnover_diff(meta_clean$SampleType)
perm <- replicate(NPERM, {
  t2 <- meta_clean$SampleType
  for (s in unique(meta_clean$Site)) {
    k <- meta_clean$Site == s
    t2[k] <- sample(t2[k])
  }
  turnover_diff(t2)
})
p_val <- (sum(perm >= obs) + 1) / (NPERM + 1)

cat("\n## River-sea turnover\n")
cat("water median:  ", median(D[cross & outer(type == "Water",   type == "Water")]), "\n")
cat("plastic median:", median(D[cross & outer(type == "Plastic", type == "Plastic")]), "\n")
cat("difference:    ", obs, "\n")
cat("permutation p: ", p_val, "\n")

# -- 7. River or sea x sample type PERMANOVA -----------------------------------
# Independent check of the interaction reported from
# analysis/01_diversity_permanova.py, with all samples permuted freely.
set.seed(42)
meta_clean$RiverSea <- ifelse(meta_clean$Site %in% c("R1", "R2"), "river", "sea")
cat("\n## River or sea x sample type PERMANOVA\n")
print(adonis2(bray ~ RiverSea * SampleType, data = meta_clean,
              by = "terms", permutations = NPERM))
