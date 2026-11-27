# Dispersion and river-sea turnover
#
# Homogeneity of multivariate dispersion (betadisper / permutest) for the
# Bray-Curtis matrix used in the PERMANOVA, and a within-site permutation test
# of the difference in river-sea turnover between plastic-bound and planktonic
# communities, and a check of the river or sea x sample type PERMANOVA.
#
# The count matrix, site relabelling and Bray-Curtis matrix are built the same
# way as for the PCoA (Bray-Curtis on species counts).
#
# Run from the repository root:
#   Rscript analysis/02_beta_diversity/dispersion_turnover.R

source("analysis/paths.R")
OUT <- step_output("02_beta_diversity")

library(vegan)
library(dplyr)

set.seed(42)
NPERM <- 9999

# -- 1. Load data and build the Bray-Curtis matrix -----------------------------
sp_tab <- seq_labels(read.csv(
  file.path(DATA, "counts.csv"),
  header = TRUE, check.names = FALSE))

# Library columns are named site + sample type + replicate, e.g. R3P2.
libs <- grep("^(R2|R3|S1|Sc)[PW][0-9]$", colnames(sp_tab), value = TRUE)
stopifnot(length(libs) == 24)
sp_counts <- sp_tab[, libs]
rownames(sp_counts) <- sp_tab$species
sp_counts_clean <- sp_counts[rowSums(sp_counts) > 0, ]

# Site and sample type come from the library names. Sequencing labels R2 and
# R3 are manuscript R1 and R2.
meta_clean <- data.frame(
  SampleName = libs,
  Site       = dplyr::recode(substr(libs, 1, 2), "R2" = "R1", "R3" = "R2"),
  SampleType = ifelse(substr(libs, 3, 3) == "P", "Plastic", "Water"),
  stringsAsFactors = FALSE)

bray <- vegdist(t(sp_counts_clean), method = "bray")
stopifnot(all(labels(bray) == meta_clean$SampleName))

site <- factor(meta_clean$Site, levels = c("R1", "R2", "S1", "Sc"))
type <- factor(meta_clean$SampleType)
rs   <- factor(ifelse(site %in% c("R1", "R2"), "river", "sea"))
print(table(site, type))
print(adonis2(bray ~ site + type, by = "terms", permutations = NPERM))

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
# analysis/02_beta_diversity/script.py, with all samples permuted freely.
set.seed(42)
meta_clean$RiverSea <- ifelse(meta_clean$Site %in% c("R1", "R2"), "river", "sea")
cat("\n## River or sea x sample type PERMANOVA\n")
print(adonis2(bray ~ RiverSea * SampleType, data = meta_clean,
              by = "terms", permutations = NPERM))
