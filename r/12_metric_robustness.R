# Sensitivity of the river-sea turnover result to the distance measure
#
# Repeats the river or sea x sample type PERMANOVA (adonis2, terms added
# sequentially, 9,999 permutations) and the within-site permutation test of
# river-sea turnover (r/08_dispersion_turnover.R) under six distance measures.
# Also reports, for each fraction, the ratio of the median river-sea distance to
# the median distance between sites within the same salinity regime.
# Reported in Results 3.2 and Supplementary-1.
#
# Run from the repository root:
#   Rscript r/12_metric_robustness.R

source("r/00_paths.R")
suppressMessages(library(vegan))
set.seed(42)
NPERM <- 9999

x <- read.csv(file.path(DATA, "counts.csv"), check.names = FALSE)
tax <- intersect(c("total", "kingdom", "phylum", "class", "order", "family",
                   "genus", "tax", "superkingdom", "species"), colnames(x))
m <- as.matrix(x[, setdiff(colnames(x), tax)]); m[is.na(m)] <- 0
X <- t(m); X <- X[, colSums(X) > 0]

# sequencing labels: R2* = R1, R3* = R2
s2   <- substr(rownames(X), 1, 2)
site <- ifelse(s2 == "R2", "R1", ifelse(s2 == "R3", "R2", s2))
type <- ifelse(substr(rownames(X), 3, 3) == "P", "Plastic", "Water")
rs   <- ifelse(site %in% c("R1", "R2"), "river", "sea")
md   <- data.frame(site, type, rs)
cat("Median species richness: plastic", median(rowSums(X[type == "Plastic", ] > 0)),
    "| water", median(rowSums(X[type == "Water", ] > 0)), "\n\n")

cross   <- outer(rs == "river", rs == "sea")
samereg <- outer(rs, rs, "==") & !outer(site, site, "==")
same_t  <- function(tp, t) outer(tp == t, tp == t)
turn    <- function(D, tp) median(D[cross & same_t(tp, "Water")]) - median(D[cross & same_t(tp, "Plastic")])
clr     <- function(M) { L <- log(M + 1); L - rowMeans(L) }

dists <- list(
  "Bray-Curtis, read counts"         = function() vegdist(X, "bray"),
  "Bray-Curtis, rarefied"            = function() vegdist(rrarefy(X, min(rowSums(X))), "bray"),
  "Bray-Curtis, relative abundance"  = function() vegdist(X / rowSums(X), "bray"),
  "Jaccard, presence-absence"        = function() vegdist(X, "jaccard", binary = TRUE),
  "Hellinger"                        = function() dist(decostand(X, "hellinger")),
  "Aitchison, CLR with pseudocount 1" = function() dist(clr(X)))

out <- do.call(rbind, lapply(names(dists), function(nm) {
  d <- dists[[nm]](); D <- as.matrix(d)
  a <- adonis2(d ~ rs * type, data = md, by = "terms", permutations = NPERM)
  obs  <- turn(D, type)
  perm <- replicate(NPERM, { t2 <- type
    for (s in unique(site)) { k <- site == s; t2[k] <- sample(t2[k]) }
    turn(D, t2) })
  med <- function(mask, t) median(D[mask & same_t(type, t)])
  data.frame(metric = nm,
    water_river_sea   = signif(med(cross, "Water"), 4),
    plastic_river_sea = signif(med(cross, "Plastic"), 4),
    turnover_p        = (sum(perm >= obs) + 1) / (NPERM + 1),
    water_ratio       = round(med(cross, "Water")   / med(samereg, "Water"), 2),
    plastic_ratio     = round(med(cross, "Plastic") / med(samereg, "Plastic"), 2),
    interaction_R2    = round(a$R2[3], 3),
    interaction_F     = round(a$F[3], 2),
    interaction_p     = a$`Pr(>F)`[3])
}))
options(width = 200)
print(out, row.names = FALSE)
write.csv(out, file.path(OUT, "metric_robustness.csv"), row.names = FALSE)
