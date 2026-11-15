#!/usr/bin/env Rscript
# =============================================================================
# Figure 5C: CTX-M group 1 amplicon tree, pruned to a panel.
#
#   Rscript phylogeny/ctxm/06_plot_figure5c.R
#
# Inputs, all in the repository:
#   output/ctxm/ctxm_group1.treefile     the deposited IQ-TREE tree
#   output/ctxm/haplotype_collapse.tsv   members of each collapsed reference tip
#   phylogeny/ctxm/figure5c_tips.txt     which tips the panel draws
#
# Steps, in this order:
#   1. Root on the group 9 outgroup, as in the Methods.
#   2. Collapse zero-length branches in the COMPLETE tree. Doing this before
#      pruning means every node left in the panel is a clade of the full tree.
#   3. Prune to the tips in figure5c_tips.txt. The outgroup is not drawn.
#   4. Support at each panel node is the SH-aLRT/UFBoot of the same clade in the
#      complete tree, looked up by its tip set, so pruning cannot misplace it.
#   5. The two branches from the root are drawn at a fixed length and marked
#      with a break, because the root position is set by the undrawn outgroup.
#
# Study sequences:
#   filled star, red box   the blaCTX-M-177-like sequence; R1W stands for the
#                          five identical consensus sequences
#   open star, "S1W"       blaCTX-M-36, identical to the S1W sequence
# =============================================================================
suppressPackageStartupMessages({
  library(ggtree); library(ggplot2); library(ape)
})

TREE    <- Sys.getenv("CTXM_TREE",   "output/ctxm/ctxm_group1.treefile")
TIPS    <- Sys.getenv("CTXM_TIPS",   "phylogeny/ctxm/figure5c_tips.txt")
OUTDIR  <- Sys.getenv("CTXM_OUTDIR", "figures")
FIG_W   <- as.numeric(Sys.getenv("CTXM_W", "6.5"))
FIG_H   <- as.numeric(Sys.getenv("CTXM_H", "2.8"))
HAPNAME <- Sys.getenv("CTXM_HAPNAME", "blaCTX-M-177-like Taninim Haplotype")
SCALE   <- Sys.getenv("CTXM_SCALE", "yes")      # "no" hides the scale bar

source("phylogeny/ctxm/plot_helpers.R")
OGPAT <- "^blaCTX-M-(14|27|65)(_and_[0-9]+_more)?$"

# --- 1-2. full tree: root, collapse --------------------------------------------
full <- ape::read.tree(TREE)
og   <- grep(OGPAT, full$tip.label, value = TRUE)
full <- ape::root(full, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)
full <- collapse_zero(full)

# --- 3. prune ------------------------------------------------------------------
keep <- trimws(readLines(TIPS, warn = FALSE))
keep <- keep[nzchar(keep) & !startsWith(keep, "#")]
miss <- setdiff(keep, full$tip.label)
if (length(miss)) stop("tips not in the tree: ", paste(miss, collapse = ", "))
tr <- ape::keep.tip(full, keep)

# --- 4. support by clade -------------------------------------------------------
n   <- length(tr$tip.label)
nf  <- length(full$tip.label)
fmt <- function(s) {
  p <- suppressWarnings(as.numeric(strsplit(s, "/")[[1]]))
  if (length(p) != 2 || anyNA(p)) return(NA_character_)
  sprintf("%.0f/%.0f", p[1], p[2])
}
sup <- vapply((n + 1):(n + tr$Nnode), function(k) {
  if (k == n + 1) return(NA_character_)            # the root: no clade to support
  tips <- ape::extract.clade(tr, k)$tip.label
  fmt(full$node.label[ape::getMRCA(full, tips) - nf])
}, character(1))
supdf <- data.frame(node = (n + 1):(n + tr$Nnode), sup = sup)

# --- 5. fixed-length root branches ---------------------------------------------
ROOT_LEN <- 0.009
root_edges <- which(tr$edge[, 1] == n + 1)
true_root <- tr$edge.length[root_edges]
tr$edge.length[root_edges] <- ROOT_LEN
message("root branches drawn at ", ROOT_LEN, "; true lengths ",
        paste(sprintf("%.4f", true_root), collapse = ", "))

# --- labels --------------------------------------------------------------------
ct <- utils::read.delim(sub("[^/]+$", "haplotype_collapse.tsv", TREE),
                        stringsAsFactors = FALSE)
members <- setNames(strsplit(ct$members, ","), ct$tip_label)
num <- function(x) as.integer(sub("^blaCTX-M-", "", x))

tip_label <- function(lab) {
  if (lab == "R1W") return(HAPNAME)
  m <- members[[lab]]
  if (is.null(m) || length(m) == 1) return(lab)
  m <- m[order(num(m))]
  if (length(m) == 2) return(sprintf("blaCTX-M-%d / -%d", num(m[1]), num(m[2])))
  sprintf("%s group (n = %d)", sub("_and_.*$", "", lab), length(m))
}
meta <- data.frame(label = tr$tip.label,
                   show  = vapply(tr$tip.label, tip_label, character(1)),
                   star  = ifelse(tr$tip.label == "R1W", "filled",
                           ifelse(tr$tip.label == "blaCTX-M-36", "open", NA)),
                   stringsAsFactors = FALSE)
meta$show <- ifelse(is.na(meta$star), meta$show, paste0("     ", meta$show))

NOTES <- data.frame(
  label = c("R1W", "blaCTX-M-36"),
  text  = c("R1W, R1P,\nR2W, S1W,\nScW", "S1W"),
  boxed = c(TRUE, FALSE))

# --- draw ----------------------------------------------------------------------
TIPSIZE <- 2.5                                        # mm
p <- ggtree(ape::ladderize(tr, right = TRUE), linewidth = 0.35) %<+% meta %<+% supdf
d <- p$data

# Label widths are MEASURED with the same font and size, then converted to data
# units, so the arrows start where each label ends whatever the label text.
xmax    <- max(d$x)
TREE_FR <- 0.42                                       # share of width for the tree
xr      <- xmax / TREE_FR
lab_off <- xmax * 0.03
MARGIN_IN <- 8 / 72                                   # plot.margin, left + right
data_per_in <- (xr * 1.04) / (FIG_W - MARGIN_IN)
grDevices::pdf(NULL)
text_in <- function(txt) grid::convertWidth(grid::grobWidth(grid::textGrob(txt,
  gp = grid::gpar(fontsize = TIPSIZE * ggplot2::.pt, fontface = "bold.italic"))),
  "inches", valueOnly = TRUE)
tipd    <- d[d$isTip, ]
tipd$lab_end <- tipd$x + lab_off +
  vapply(tipd$show, text_in, numeric(1)) * data_per_in
star_x  <- text_in("  ") * data_per_in              # star sits in the leading gap
invisible(grDevices::dev.off())

ann <- merge(tipd[, c("label", "y", "lab_end")], NOTES, by = "label")
arrow_to <- max(ann$lab_end) + xr * 0.08
ann$x0 <- ann$lab_end + xr * 0.012
ann$x1 <- arrow_to

p <- p +
  geom_tiplab(aes(label = show), fontface = "bold.italic", size = TIPSIZE,
              offset = lab_off) +
  # five-point stars from DejaVu Sans; ggplot has no star point shape
  geom_text2(aes(subset = isTip & !is.na(star), x = x + lab_off + star_x,
                 label = ifelse(star == "filled", "\u2605", "\u2606")),
             family = "DejaVu Sans", size = 3.2, hjust = 0.5) +
  geom_text2(aes(subset = !isTip & !is.na(sup), label = sup),
             fontface = "bold", size = 2.2, hjust = 1.08, vjust = -0.35) +
  geom_segment(data = ann, aes(x = x0, xend = x1, y = y, yend = y),
               inherit.aes = FALSE, colour = "red", linewidth = 0.4,
               arrow = arrow(length = unit(1.4, "mm"), type = "open")) +
  geom_label(data = ann[ann$boxed, ], aes(x = x1 + xr * 0.01, y = y, label = text),
             inherit.aes = FALSE, hjust = 0, size = 2.3, fontface = "bold",
             colour = "black", border.colour = "red", linewidth = 0.5,
             label.r = unit(0, "pt"), label.padding = unit(0.3, "lines"),
             lineheight = 0.95) +
  geom_text(data = ann[!ann$boxed, ], aes(x = x1 + xr * 0.01, y = y, label = text),
            inherit.aes = FALSE, hjust = 0, size = 2.3, fontface = "bold") +
  annotate("text", x = -xr * 0.03, y = max(d$y) + 0.5, label = "C",
           fontface = "bold", size = 4, hjust = 0)

# break marks on the two root branches
rb <- d[d$parent == n + 1 & d$node != n + 1, ]
rx <- d$x[d$node == n + 1] + ROOT_LEN * 0.3
bw <- ROOT_LEN * 0.07
brk <- do.call(rbind, lapply(rb$y, function(y)
  data.frame(x = c(rx - bw, rx + bw) - bw, xend = c(rx, rx + 2 * bw) - bw,
             y = y - 0.3, yend = y + 0.3)))
p <- p +
  geom_rect(data = data.frame(y = rb$y), aes(xmin = rx - 1.6 * bw, xmax = rx + 1.2 * bw,
            ymin = y - 0.12, ymax = y + 0.12), inherit.aes = FALSE,
            fill = "white", colour = NA) +
  geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.35)

if (SCALE == "yes")
  p <- p + geom_treescale(x = xmax * 0.6, y = -0.4, width = 0.005,
                          fontsize = 2.1, offset = 0.3, linesize = 0.35)

p <- p + theme_tree() +
  coord_cartesian(xlim = c(-xr * 0.04, xr), ylim = c(-0.8, max(d$y) + 0.8), expand = FALSE,
                  clip = "off") +
  theme(plot.margin = margin(4, 4, 4, 4),
        plot.background = element_rect(fill = "white", colour = NA))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
pdf_out <- file.path(OUTDIR, "Figure5C_CTXM.pdf")
png_out <- file.path(OUTDIR, "Figure5C_CTXM.png")
ggsave(pdf_out, p, width = FIG_W, height = FIG_H, device = cairo_pdf)  # cairo: star glyphs
ggsave(png_out, p, width = FIG_W, height = FIG_H, dpi = 600)
cat(pdf_out, " / ", png_out, " written\n", sep = "")
cat("support (SH-aLRT/UFBoot) by clade of the complete tree:\n")
print(supdf[!is.na(supdf$sup), ], row.names = FALSE)
