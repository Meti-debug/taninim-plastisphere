#!/usr/bin/env Rscript
# =============================================================================
# Figure 5 panel C: CTX-M group 1 tree, drawn for a journal panel.
#
#   Rscript phylogeny/ctxm/06_plot_figure5c.R
#
# Differences from the general plotting scripts, all of them deliberate:
#
#   * Sequences that are IDENTICAL over the amplified window are drawn as ONE
#     tip. R1W, R1P, R2W and ScW sit at IQ-TREE's minimum branch length from
#     each other, so drawing four tips would imply divergence that is not in
#     the data. S1W is NOT collapsed with them: it is one substitution away.
#   * The group 9 outgroup is dropped. Its branch is about 0.5 where the whole
#     ingroup spans 0.01, so including it flattens the ingroup into a line at
#     panel size. The caption states the rooting; the full tree shows it.
#   * No support values. The one well-supported node, group 1 versus group 9 at
#     100/100, spans all 68 group 1 tips and is not a clade this panel draws, so
#     no number here could stand for it.
#   * Tree lines are black. Colour is used only for the study annotations.
# =============================================================================
suppressPackageStartupMessages({
  library(ggtree); library(ggplot2); library(ape); library(dplyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

TREE    <- Sys.getenv("CTXM_TREE",   "output/ctxm/ctxm_group1.treefile")
OUTDIR  <- Sys.getenv("CTXM_OUTDIR", "figures")
KEEP_N  <- as.integer(Sys.getenv("CTXM_KEEP", "8"))
FIG_W   <- as.numeric(Sys.getenv("CTXM_W", "3.2"))
FIG_H   <- as.numeric(Sys.getenv("CTXM_H", "2.6"))
HAPNAME <- Sys.getenv("CTXM_HAPNAME", "blaCTX-M-177-like Taninim haplotype")
VARNAME <- Sys.getenv("CTXM_VARNAME", "Taninim variant")
NOTEX   <- as.numeric(Sys.getenv("CTXM_NOTEX", "3.0"))  # red box column, as a
# multiple of the tree's x extent. Tip labels are drawn in device space and run
# well past the tree, so this has to clear the longest label; raise it if the
# boxes sit on the text, lower it if they drift off the panel.

if (!file.exists(TREE)) stop("tree not found: ", TREE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

tr <- ape::read.tree(TREE)
OGPAT <- "^blaCTX-M-(14|27|65)(_and_[0-9]+_more)?$"

# Collapsed reference tips are written by 01 as blaCTX-M-15_and_50_more. Render
# them the way the published figure does: two alleles as "blaCTX-M-12 / -96",
# more than two as "blaCTX-M-15 group (n = 51)". haplotype_collapse.tsv, also
# written by 01, holds the members, so pairs can name both alleles exactly.
COLLAPSE_TSV <- file.path(dirname(TREE), "haplotype_collapse.tsv")
members_of <- list()
if (file.exists(COLLAPSE_TSV)) {
  ct <- utils::read.delim(COLLAPSE_TSV, stringsAsFactors = FALSE)
  members_of <- setNames(strsplit(ct$members, ","), ct$tip_label)
  message("collapse table read: ", nrow(ct), " reference tips")
} else {
  message("no haplotype_collapse.tsv; pair labels will fall back to counts")
}

group_label <- function(lab) {
  m <- regmatches(lab, regexec("^(.*)_and_([0-9]+)_more$", lab))[[1]]
  if (!length(m)) return(sub("_and_(\\d+)_more", "", lab))
  head_allele <- m[2]; n <- as.integer(m[3]) + 1
  mem <- members_of[[lab]]
  if (n == 2 && length(mem) == 2) {
    other <- sub("^blaCTX-M", "", setdiff(mem, head_allele)[1])
    return(paste0(head_allele, " /", other))
  }
  paste0(head_allele, " group (n = ", n, ")")
}

TIPFILE <- Sys.getenv("CTXM_STUDY_TIPS",
                      file.path(dirname(TREE), "study_tips.txt"))
study_tips <- trimws(readLines(TIPFILE, warn = FALSE))
study_tips <- study_tips[nzchar(study_tips)]

# --- which study tips are identical to each other over the window? -----------
# cophenetic distance of 0 means no substitution separates them.
d <- ape::cophenetic.phylo(tr)
site_tips <- intersect(c("R1W","R1P","R2W","R2P","S1W","S1P","ScW","ScP"), study_tips)
groups <- list()
for (t in site_tips) {
  placed <- FALSE
  for (g in seq_along(groups)) {
    if (d[t, groups[[g]][1]] < 1e-5) { groups[[g]] <- c(groups[[g]], t); placed <- TRUE; break }
  }
  if (!placed) groups[[length(groups) + 1]] <- t
}
for (g in groups)
  message("identical over the window: ", paste(g, collapse = ", "))

groups <- groups[order(vapply(groups, length, integer(1)), decreasing = TRUE)]
main_group <- groups[[1]]          # the shared haplotype
reps <- vapply(groups, `[`, character(1), 1)      # one tip drawn per group
collapsed_away <- setdiff(site_tips, reps)

# --- which tips to draw ------------------------------------------------------
og_tips <- grep(OGPAT, tr$tip.label, value = TRUE)
other_study <- setdiff(study_tips, site_tips)      # e.g. blaCTX-M-36_S1W
anchor <- c(reps, other_study)
cand <- setdiff(tr$tip.label, c(study_tips, og_tips))
nearest <- sort(apply(d[cand, anchor, drop = FALSE], 1, min))
keep <- unique(c(anchor, names(head(nearest, KEEP_N))))
sub <- ape::keep.tip(tr, keep)
message("drawing ", length(keep), " tips (", length(collapsed_away),
        " collapsed, ", length(og_tips), " outgroup dropped)")

# --- labels and the red sample annotations -----------------------------------
SITELAB <- c(R1W = "R1W", R1P = "R1P", R2W = "R2W", R2P = "R2P",
             S1W = "S1W", S1P = "S1P", ScW = "ScW", ScP = "ScP")

tip_show <- character(length(sub$tip.label))
tip_note <- rep(NA_character_, length(sub$tip.label))
is_study <- logical(length(sub$tip.label))

for (i in seq_along(sub$tip.label)) {
  lab <- sub$tip.label[i]
  g <- Filter(function(x) lab %in% x, groups)
  if (length(g)) {                                   # a collapsed site group
    members <- g[[1]]
    is_study[i] <- TRUE
    # only the shared haplotype carries the haplotype name; a group that is a
    # substitution away from it is a different sequence and is named as one
    tip_show[i] <- if (identical(members, main_group)) HAPNAME else VARNAME
    tip_note[i] <- paste(SITELAB[members], collapse = ", ")
  } else if (lab %in% other_study) {                 # e.g. blaCTX-M-36_S1W
    is_study[i] <- TRUE
    tip_show[i] <- sub("_[^_]+$", "", lab)
    tip_note[i] <- sub("^.*_", "", lab)
  } else {
    tip_show[i] <- group_label(lab)
  }
}

meta <- data.frame(label = sub$tip.label, show = tip_show,
                   note = tip_note, study = is_study,
                   stringsAsFactors = FALSE)

p <- ggtree(sub, size = 0.5, colour = "black") %<+% meta

# Place the red boxes by computing how much data-space the tip labels occupy.
#
# Guessing an offset cannot work: tip labels are drawn at a fixed physical size,
# so widening the panel shrinks the tree in data units and the text takes up
# MORE data-space. Raising the offset and the expansion together just lets the
# label chase the box. Instead, convert the font size to inches, divide by the
# panel width to get a fraction, and solve for the expansion that fits both the
# longest label and the longest note.
TIPSIZE <- 2.6                       # ggplot size unit = mm
char_in <- 0.60 * TIPSIZE / 25.4     # a proportional font averages ~0.6 of its height
nmax    <- max(nchar(meta$show))
notemax <- max(nchar(meta$note), na.rm = TRUE)

xmax    <- max(p$data$x, na.rm = TRUE)
lab_off <- xmax * 0.04
K       <- nmax + notemax + 4        # characters to fit right of the tree
a       <- K * char_in / FIG_W
if (a >= 0.95)
  stop("labels need more width than the panel has; raise CTXM_W or cut CTXM_KEEP")
extra   <- (lab_off + a * xmax) / (1 - a)
char_x  <- char_in / FIG_W * (xmax + extra)
note_x  <- xmax + lab_off + (nmax + 1) * char_x
message(sprintf("panel: %d tips, longest label %d chars, note column at %.1f%% of width",
                length(keep), nmax, 100 * note_x / (xmax + extra)))

p <- p +
  geom_tiplab(aes(label = show, fontface = ifelse(study, 2, 3)),
              size = TIPSIZE, offset = lab_off, align = TRUE,
              linetype = "dotted", linesize = 0.2, colour = "black",
              show.legend = FALSE) +
  geom_tippoint(size = 1.1, colour = "black") +
  geom_label(data = function(d) d[!is.na(d$note) & d$isTip, ],
             aes(x = note_x, label = note),
             size = 2.4, colour = "red", fill = "white",
             label.size = 0.3, label.padding = unit(0.10, "lines"),
             hjust = 0, show.legend = FALSE) +
  geom_treescale(width = 0.002, fontsize = 2.3, offset = 0.4) +
  theme_tree() +
  ggplot2::expand_limits(x = xmax + extra) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(4, 6, 4, 4))

pdf_out <- file.path(OUTDIR, "Figure5C_CTXM.pdf")
png_out <- file.path(OUTDIR, "Figure5C_CTXM.png")
ggsave(pdf_out, p, width = FIG_W, height = FIG_H, limitsize = FALSE)
ggsave(png_out, p, width = FIG_W, height = FIG_H, dpi = 600, limitsize = FALSE)
cat(pdf_out, " / ", png_out, " written (", length(keep), " tips, ",
    FIG_W, " x ", FIG_H, " in)\n", sep = "")
