#!/usr/bin/env Rscript
# =============================================================================
# Figure S6, panels A and B.
#
#   Rscript analysis/phylogeny/ctxm/04_plot_figure_s6.R
#
# A  the same tree as Figure 5C (panel_tree() in plot_helpers.R)
# B  the variable sites among the drawn sequences, read from the alignment the
#    tree was inferred from (output/ctxm/tree_input.fasta). Rows follow the tree.
#    The blaCTX-M-177-like variant is shown in full; other rows show a base only
#    where it differs from the variant. "diff." is the number of differences.
#    Columns are CDS positions (alignment column 1 = CDS position 28).
# =============================================================================
suppressPackageStartupMessages({
  library(ggtree); library(ggplot2); library(ape)
})

TREE    <- Sys.getenv("CTXM_TREE",   "output/ctxm/ctxm_group1.treefile")
TIPS    <- Sys.getenv("CTXM_TIPS",   "analysis/phylogeny/ctxm/figure5c_tips.txt")
ALN     <- Sys.getenv("CTXM_ALN",    "output/ctxm/tree_input.fasta")
OUTDIR  <- Sys.getenv("CTXM_OUTDIR", "output/figures")
FIG_W   <- as.numeric(Sys.getenv("CTXM_W", "8"))
FIG_H   <- as.numeric(Sys.getenv("CTXM_H", "4.6"))
WIN_LO  <- 28                                   # CDS position of alignment column 1
FOCAL   <- 46                                   # boxed position

source("analysis/phylogeny/ctxm/plot_helpers.R")
pt <- panel_tree(TREE, TIPS, ROOT_LEN = 0.009, HAPNAME = "blaCTX-M-177-like")
tr <- pt$tree; n <- pt$n

# --- alignment, variable sites --------------------------------------------------
read_fasta <- function(f) {
  l <- readLines(f); h <- grep("^>", l)
  s <- vapply(seq_along(h), function(i) {
    e <- if (i < length(h)) h[i + 1] - 1 else length(l)
    toupper(paste(l[(h[i] + 1):e], collapse = ""))
  }, character(1))
  setNames(s, sub("^>(\\S+).*", "\\1", l[h]))
}
aln <- read_fasta(ALN)
miss <- setdiff(tr$tip.label, names(aln))
if (length(miss)) stop("tips missing from the alignment: ", paste(miss, collapse = ", "))
M <- do.call(rbind, strsplit(aln[tr$tip.label], ""))
rownames(M) <- tr$tip.label
var <- which(apply(M, 2, function(col) length(setdiff(unique(col), c("N", "-"))) > 1))
ref <- M["R1W", ]

# --- tree panel -----------------------------------------------------------------
meta <- data.frame(label = tr$tip.label,
                   show  = vapply(tr$tip.label, pt$tip_label, character(1)),
                   star  = ifelse(tr$tip.label == "R1W", "filled",
                           ifelse(tr$tip.label == "blaCTX-M-36", "open", NA)),
                   stringsAsFactors = FALSE)
meta$show[meta$label == "blaCTX-M-36"] <- "blaCTX-M-36 (S1 only)"
meta$show <- ifelse(is.na(meta$star), meta$show, paste0("     ", meta$show))

TIPSIZE <- 2.6
p <- ggtree(ape::ladderize(tr, right = TRUE), linewidth = 0.35) %<+% meta %<+% pt$support
d <- p$data
xmax <- max(d$x)

# physical layout, in inches across the figure width
TREE_IN <- 1.9
dpi_x   <- xmax / TREE_IN                      # data units per inch
grDevices::pdf(NULL)
text_in <- function(txt, size = TIPSIZE, face = "bold.italic")
  grid::convertWidth(grid::grobWidth(grid::textGrob(txt,
    gp = grid::gpar(fontsize = size * ggplot2::.pt, fontface = face))),
    "inches", valueOnly = TRUE)
lab_off <- xmax * 0.03
lab_end <- max(d$x[d$isTip]) + lab_off +
  max(vapply(meta$show, text_in, numeric(1))) * dpi_x
star_x  <- text_in("  ") * dpi_x
invisible(grDevices::dev.off())

# --- matrix panel ---------------------------------------------------------------
CELL_IN <- 0.2
cw  <- CELL_IN * dpi_x
x0  <- lab_end + 0.25 * dpi_x                  # left edge of the first column
tip_y <- setNames(d$y[d$isTip], d$label[d$isTip])
cols <- data.frame(k = seq_along(var), pos = WIN_LO + var - 1)
cols$x <- x0 + (cols$k - 0.5) * cw

cells <- do.call(rbind, lapply(rownames(M), function(t) {
  b <- M[t, var]
  is_ref <- t == "R1W"
  show <- is_ref | b != ref[var]
  data.frame(tip = t, k = seq_along(var), base = b, show = show,
             fillv = if (is_ref) rep("ref", length(b)) else b, y = tip_y[[t]])
}))
cells$x <- x0 + (cells$k - 0.5) * cw
cells <- cells[cells$show, ]

ndiff <- data.frame(tip = rownames(M), y = tip_y[rownames(M)],
                    nd = vapply(rownames(M), function(t) sum(M[t, var] != ref[var]), numeric(1)))
diff_x <- x0 + (length(var) + 0.9) * cw

BASECOL <- c(A = "#00A651", C = "#0072BC", G = "#F7941D", T = "#ED1C24", ref = "#FFD54F")
ymax <- max(d$y)
k46  <- which(cols$pos == FOCAL)

p <- p +
  geom_tiplab(aes(label = show), fontface = "bold.italic", size = TIPSIZE,
              offset = lab_off) +
  geom_text2(aes(subset = isTip & !is.na(star), x = x + lab_off + star_x,
                 label = ifelse(star == "filled", "\u2605", "\u2606")),
             family = "DejaVu Sans", size = 3.4, hjust = 0.5) +
  geom_text2(aes(subset = !isTip & !is.na(sup), label = sup),
             fontface = "bold", size = 2.3, hjust = 1.08, vjust = -0.35) +
  # matrix
  geom_tile(data = cells, aes(x = x, y = y, fill = fillv), inherit.aes = FALSE,
            width = cw * 0.92, height = 0.8) +
  geom_text(data = cells, aes(x = x, y = y, label = base), inherit.aes = FALSE,
            size = 2.5, fontface = "bold") +
  scale_fill_manual(values = BASECOL, guide = "none") +
  annotate("rect", xmin = cols$x[k46] - cw / 2, xmax = cols$x[k46] + cw / 2,
           ymin = 0.45, ymax = ymax + 0.55, fill = NA, colour = "black", linewidth = 0.5) +
  geom_text(data = cols, aes(x = x, y = ymax + 0.9, label = pos), inherit.aes = FALSE,
            angle = 90, hjust = 0, size = 2.4, fontface = "bold") +
  annotate("text", x = diff_x, y = ymax + 0.9, label = "diff.", angle = 90, hjust = 0,
           size = 2.4, fontface = "bold") +
  geom_text(data = ndiff, aes(x = diff_x, y = y, label = nd), inherit.aes = FALSE,
            size = 2.5, fontface = "bold") +
  annotate("segment", x = cols$x[k46], xend = cols$x[k46], y = ymax + 3.3, yend = ymax + 2.5,
           arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.4) +
  annotate("text", x = cols$x[k46] + 0.3 * cw, y = ymax + 3.4,
           label = "nt 46 A>G (Thr16Ala)", hjust = 0, size = 2.6, fontface = "bold")

# break marks on the two root branches
ROOT_LEN <- 0.009
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

# key: stars and scale, top left; panel letters
p <- p +
  annotate("text", x = 0, y = ymax + 3.4, label = "\u2605", family = "DejaVu Sans",
           size = 3.6, hjust = 0) +
  annotate("text", x = 0.22 * dpi_x, y = ymax + 3.4, label = "this study: detected at all sites",
           size = 2.5, hjust = 0, fontface = "bold") +
  annotate("text", x = 0, y = ymax + 2.6, label = "\u2606", family = "DejaVu Sans",
           size = 3.6, hjust = 0) +
  annotate("text", x = 0.22 * dpi_x, y = ymax + 2.6, label = "this study: detected at S1 only",
           size = 2.5, hjust = 0, fontface = "bold") +
  annotate("segment", x = 0, xend = 0.0027, y = ymax + 1.8, yend = ymax + 1.8, linewidth = 0.5) +
  annotate("text", x = 0.0027 + 0.05 * dpi_x, y = ymax + 1.8, hjust = 0, size = 2.5,
           fontface = "bold", label = "0.0027 subs/site (\u22481 nt)") +
  annotate("text", x = -0.1 * dpi_x, y = ymax + 0.6, label = "A", fontface = "bold",
           size = 5, hjust = 0) +
  annotate("text", x = x0 - 0.15 * dpi_x, y = ymax + 3.4, label = "B", fontface = "bold",
           size = 5, hjust = 1)

total_in <- (diff_x + 0.3 * dpi_x) / dpi_x
p <- p + theme_tree() +
  coord_cartesian(xlim = c(-0.15 * dpi_x, diff_x + 0.3 * dpi_x),
                  ylim = c(0.3, ymax + 3.9), expand = FALSE, clip = "off") +
  theme(plot.margin = margin(4, 4, 4, 4),
        plot.background = element_rect(fill = "white", colour = NA))

W <- total_in + 0.25
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
pdf_out <- file.path(OUTDIR, "FigureS6_AB.pdf")
png_out <- file.path(OUTDIR, "FigureS6_AB.png")
ggsave(pdf_out, p, width = W, height = FIG_H, device = cairo_pdf)
ggsave(png_out, p, width = W, height = FIG_H, dpi = 600)
cat(pdf_out, " / ", png_out, " written (", round(W, 2), " x ", FIG_H, " in)\n", sep = "")
cat("variable CDS positions:", paste(cols$pos, collapse = " "), "\n")
