#!/usr/bin/env Rscript
# =============================================================================
# CTX-M phylogeny with side panels
#   panel 1  tree; tips = gene name + sample code
#   panel 2  sampling location and matrix
#   panel 3  horizontal abundance bars
#
# NOTE ON THE CODES: over positions called in all seven samples the sequences
# are identical. The codes identify samples, not distinct alleles, and the
# branching among them reflects which positions were masked (N) in each
# consensus rather than real sequence divergence. Nothing in the topology
# among the study tips supports a claim about transport along the transect.
# See README.md.
#
# EXPECTED INPUT LABEL FORMAT
#   samples     CTXM15_<site>_<matrix>_bc<NN>
#   references  blaCTX-M-<n>___g1
#
#   export CTXM_TREE=output/ctxm/tree2.treefile
#   Rscript phylogeny/ctxm/03_plot_tree_bars.R
#
# Read counts are read from phylogeny/ctxm/ctxm_read_counts.csv so that the
# figure and the manuscript tables cannot drift apart. Regenerate that file
# from the detection step rather than editing numbers here.
# =============================================================================
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2); library(dplyr); library(readr)
})

TREE   <- Sys.getenv("CTXM_TREE",   "output/ctxm/tree2.treefile")
COUNTS <- Sys.getenv("CTXM_COUNTS", "phylogeny/ctxm/ctxm_read_counts.csv")
OUTDIR <- Sys.getenv("CTXM_OUTDIR", "figures")

GENE  <- "blaCTX-M-177-like"
VALUE <- Sys.getenv("CTXM_METRIC", "lib_pct")   # "lib_pct" or "rel_pct"

# Ingroup rooting; see the note in 02_plot_tree_gene.R and README.md.
OUTGROUP <- Sys.getenv("CTXM_OUTGROUP", "blaCTX-M-36___g1")

if (!file.exists(TREE))   stop("tree not found: ", TREE)
if (!file.exists(COUNTS)) stop("read counts not found: ", COUNTS)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

tr <- read.tree(TREE)
if (OUTGROUP %in% tr$tip.label) {
  tr <- root(tr, outgroup = OUTGROUP, resolve.root = TRUE)
  message("rooted on ", OUTGROUP)
} else {
  message("outgroup ", OUTGROUP, " not present; tree left as read (unrooted)")
}

counts <- read_csv(COUNTS, show_col_types = FALSE) %>%
  mutate(rel_pct = 100 * ctxm / sum(ctxm),
         lib_pct = 100 * ctxm / total)
stopifnot(VALUE %in% names(counts))

old <- tr$tip.label
bc  <- ifelse(grepl("^CTXM15_", old), sub(".*_(bc[0-9]+)$", "\\1", old), NA_character_)
code <- counts$code[match(bc, counts$bc)]
tr$tip.label <- ifelse(is.na(bc), sub("___g1$", "", old),
                       paste0(GENE, " [", code, "]"))
lab <- tr$tip.label

info <- data.frame(id = lab, bc = bc, stringsAsFactors = FALSE) %>%
  left_join(counts, by = "bc") %>%
  mutate(study = !is.na(bc), value = .data[[VALUE]])

p <- ggtree(tr, size = 0.6) +
  geom_tiplab(aes(fontface = ifelse(label %in% info$id[info$study],
                                    "bold", "italic")),
              size = 3.3, offset = 0.0005, align = TRUE,
              linetype = "dotted", linesize = 0.25, colour = "grey15") +
  geom_tippoint(size = 2, colour = "grey30") +
  geom_treescale(width = 0.002, fontsize = 3, offset = 0.4) +
  ggtree::hexpand(0.75)

loc  <- info %>% filter(study) %>%
        transmute(id, lab = paste0(site, " · ", matrix))
bars <- info %>% filter(study) %>% transmute(id, value, matrix)
BARMAX <- max(bars$value) * 1.35      # headroom for the value labels

PANEL <- "blaCTX-M reads (% of library)"

p2 <- facet_plot(p, panel = "Sampling location", data = loc,
                 geom = geom_text, mapping = aes(x = 0, label = lab),
                 hjust = 0, size = 3.1, colour = "grey15")

p3 <- facet_plot(p2, panel = PANEL,
                 data = bars, geom = geom_segment,
                 mapping = aes(x = 0, xend = value, y = y, yend = y,
                               colour = matrix),
                 linewidth = 5)

p4 <- facet_plot(p3, panel = PANEL,
                 data = bars, geom = geom_text,
                 mapping = aes(x = value, label = sprintf("%.3f", value)),
                 hjust = -0.15, size = 2.7, colour = "grey30")

p4 <- p4 +
  ggtree::xlim_expand(c(0, BARMAX), panel = PANEL) +
  scale_colour_manual(values = c("Water" = "#1F6FB4",
                                 "Plastic biofilm" = "#C0392B"), name = NULL) +
  theme_tree2() +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 9),
        strip.background = element_blank(),
        strip.text       = element_text(size = 9.5, face = "bold"),
        panel.spacing    = unit(0.3, "lines"),
        plot.margin      = margin(8, 24, 8, 8))

pdf_out <- file.path(OUTDIR, "Figure_CTXM_tree_bars.pdf")
png_out <- file.path(OUTDIR, "Figure_CTXM_tree_bars.png")
ggsave(pdf_out, p4, width = 12.5, height = 5.5, device = cairo_pdf)
ggsave(png_out, p4, width = 12.5, height = 5.5, dpi = 600)
cat(pdf_out, " / ", png_out, " written (metric: ", VALUE, ")\n", sep = "")
cat("\nSample codes: R1W/R1P = River 1 water/plastic, R2W = River 2 water,\n")
cat("S1W/S1P = Sea 1 water/plastic, ScW/ScP = sea control water/plastic.\n")
