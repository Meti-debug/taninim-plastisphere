#!/usr/bin/env Rscript
# =============================================================================
# CTX-M group 1 phylogeny — publication figure.
#
# Tip labels give the GENE first, with sampling details in parentheses:
#     blaCTX-M-177-like (River 1, plastic biofilm)
#     blaCTX-M-212
#
# EXPECTED INPUT LABEL FORMAT
#   samples     CTXM15_<site>_<matrix>_bc<NN>
#   references  blaCTX-M-<n>___g1
# This is the labelling of the tree used for the published figure. The tree
# written by 01_build_panel_and_tree.sh relabels tips at the alignment stage
# instead, so use 04_plot_sample_tree.R for that one. The two are not
# interchangeable; see README.md.
#
#   export CTXM_TREE=output/ctxm/tree2.treefile
#   Rscript phylogeny/ctxm/02_plot_tree_gene.R
# =============================================================================
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2)
})

TREE   <- Sys.getenv("CTXM_TREE", "output/ctxm/tree2.treefile")
OUTDIR <- Sys.getenv("CTXM_OUTDIR", "figures")

# What to call the study allele. All seven samples carry the same sequence over
# the positions called in every sample; see README.md before changing this.
STUDY_GENE <- "blaCTX-M-177-like"

# Rooting. CTX-M-36 is a GROUP-1 allele, so this is an ingroup rooting: it
# orients the drawing but does not locate the root of group 1. Rooting on a
# group-9 outgroup is the defensible choice and is what 01 does. Whichever is
# used, the Methods must say which, and must not claim both.
OUTGROUP <- Sys.getenv("CTXM_OUTGROUP", "blaCTX-M-36___g1")

if (!file.exists(TREE))
  stop("tree not found: ", TREE, "\nset CTXM_TREE to the .treefile")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

tr <- read.tree(TREE)
if (OUTGROUP %in% tr$tip.label) {
  tr <- root(tr, outgroup = OUTGROUP, resolve.root = TRUE)
  message("rooted on ", OUTGROUP)
} else {
  message("outgroup ", OUTGROUP, " not present; tree left as read (unrooted)")
}

SITE <- c(R1 = "River 1", R2 = "River 2", S1 = "Sea 1", Sc = "Sea control")

old <- tr$tip.label
new <- character(length(old)); cls <- character(length(old))

for (i in seq_along(old)) {
  v <- old[i]
  if (grepl("^CTXM15_", v)) {
    p  <- sub("^CTXM15_", "", sub("_bc[0-9]+$", "", v))
    pt <- strsplit(p, "_")[[1]]
    site    <- if (pt[1] %in% names(SITE)) SITE[[pt[1]]] else pt[1]
    plastic <- length(pt) > 1 && pt[2] == "plastic"
    new[i]  <- sprintf("%s (%s, %s)", STUDY_GENE, site,
                       if (plastic) "plastic biofilm" else "water")
    cls[i]  <- if (plastic) "This study: plastic biofilm" else "This study: water"
  } else {
    new[i] <- sub("___g1$", "", v)
    cls[i] <- if (v == "blaCTX-M-177___g1") "Closest reference (CTX-M-177)"
              else "Reference allele"
  }
}
tr$tip.label <- new
cat("tip labels:\n"); print(new)

tr  <- groupOTU(tr, split(new, cls))
pal <- c("This study: plastic biofilm"   = "#C0392B",
         "This study: water"             = "#1F6FB4",
         "Closest reference (CTX-M-177)" = "#D97D0D",
         "Reference allele"              = "#2C3E50")
shp <- c("This study: plastic biofilm" = 17, "This study: water" = 16,
         "Closest reference (CTX-M-177)" = 18, "Reference allele" = 1)

p <- ggtree(tr, aes(colour = group), size = 0.6) +
  geom_tiplab(aes(colour = group,
                  fontface = ifelse(grepl("\\(", label), "bold", "italic")),
              size = 3.4, offset = 0.0004, align = TRUE,
              linetype = "dotted", linesize = 0.25) +
  geom_tippoint(aes(colour = group, shape = group), size = 3) +
  scale_colour_manual(values = pal, breaks = names(pal), name = NULL) +
  scale_shape_manual(values  = shp, breaks = names(shp), name = NULL) +
  geom_treescale(width = 0.002, fontsize = 3.2, offset = 0.5, linesize = 0.5) +
  theme_tree() +
  theme(legend.position   = c(0.26, 0.88),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.text       = element_text(size = 9.5),
        legend.key.height = unit(0.48, "cm"),
        plot.margin       = margin(10, 130, 20, 10)) +
  ggtree::hexpand(1.15) + coord_cartesian(clip = "off")

pdf_out <- file.path(OUTDIR, "Figure_CTXM_phylogeny.pdf")
png_out <- file.path(OUTDIR, "Figure_CTXM_phylogeny.png")
ggsave(pdf_out, p, width = 10, height = 5.5, device = cairo_pdf)
ggsave(png_out, p, width = 10, height = 5.5, dpi = 600)
cat("\n", pdf_out, " / ", png_out, " written\n", sep = "")
