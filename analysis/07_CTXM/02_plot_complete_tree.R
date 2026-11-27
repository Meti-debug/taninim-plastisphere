#!/usr/bin/env Rscript
# =============================================================================
# CTX-M group 1 tree, complete: all tips, including the group 9 outgroup.
#
# EXPECTED INPUT LABEL FORMAT
#   samples     R1W, R1P, R2W, S1W, ScW, ...   (relabelled during alignment)
#   references  blaCTX-M-<n>, or blaCTX-M-<n>_and_<k>_more for collapsed tips
# This is what 01_build_panel_and_tree.sh writes.
#
# Display:
#   * zero-length branches are collapsed to polytomies, so support is shown
#     only on branches that carry at least one substitution;
#   * the group 1 stem (about 0.47) is shortened for display and marked with
#     a break and its true length.
#
# Study tips are read from study_tips.txt, written by 01_build_panel_and_tree.sh.
#
#   Rscript analysis/07_CTXM/02_plot_complete_tree.R
#
# The group 9 outgroup is drawn. DROP_OUTGROUP <- TRUE hides it.
# =============================================================================
suppressPackageStartupMessages({
  library(ggtree); library(ggplot2); library(ape); library(dplyr)
})

# base R gained %||% in 4.4; rlang exports it but dplyr does not re-export it,
# so define it here rather than depending on which of those is attached.
`%||%` <- function(a, b) if (is.null(a)) b else a

TREE   <- Sys.getenv("CTXM_TREE",   "analysis/07_CTXM/tree/ctxm_group1.treefile")
OUTDIR <- Sys.getenv("CTXM_OUTDIR", "analysis/07_CTXM/output")

DROP_OUTGROUP <- FALSE           # FALSE = show the group-9 outgroup

if (!file.exists(TREE)) stop("tree not found: ", TREE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ape::read.tree keeps the internal-node labels in $node.label, where IQ-TREE
# writes "SH-aLRT/UFBoot". Reading through treeio and then pruning loses them.
tr <- ape::read.tree(TREE)
OGPAT <- "^blaCTX-M-(14|27|65)(_and_[0-9]+_more)?$"
if (DROP_OUTGROUP) {
  drop <- grep(OGPAT, tr$tip.label, value = TRUE)
  if (length(drop)) {
    tr <- ape::drop.tip(tr, drop)
    message("outgroup tips dropped after rooting: ", paste(drop, collapse = ", "))
  }
} else {
  message("outgroup kept in the drawing: ",
          paste(grep(OGPAT, tr$tip.label, value = TRUE), collapse = ", "))
}
source("analysis/07_CTXM/plot_helpers.R")
tr <- collapse_zero(tr)
og_tips <- grep(OGPAT, tr$tip.label, value = TRUE)
STEM <- NULL
if (!DROP_OUTGROUP && length(og_tips)) {
  STEM_SUP <- stem_support(tr, og_tips)
  STEM <- shorten_stem(tr, og_tips, show_len = 0.004)
  tr <- STEM$tree
}

message("tips: ", length(tr$tip.label),
        " | internal nodes carrying a support label: ",
        sum(nzchar(tr$node.label %||% character(0))))

# Which tips are this study's is read from study_tips.txt, written by
# 01_build_panel_and_tree.sh, rather than hardcoded here. An extra haplotype
# supplied through CTXM_EXTRA is then coloured automatically instead of being
# silently drawn as a reference allele.
SITEMAP <- tibble::tribble(
  ~label, ~location,     ~matrix,
  "R1W",  "River 1",     "water",
  "R1P",  "River 1",     "plastic",
  "R2W",  "River 2",     "water",
  "R2P",  "River 2",     "plastic",
  "S1W",  "Sea 1",       "water",
  "S1P",  "Sea 1",       "plastic",
  "ScW",  "Sea control", "water",
  "ScP",  "Sea control", "plastic"
)

TIPFILE <- Sys.getenv("CTXM_STUDY_TIPS",
                      file.path(dirname(TREE), "study_tips.txt"))
study_tips <- if (file.exists(TIPFILE)) {
  trimws(readLines(TIPFILE, warn = FALSE))
} else {
  message("study_tips.txt not found at ", TIPFILE,
          "; falling back to the site codes present in the tree")
  intersect(SITEMAP$label, tr$tip.label)
}
study_tips <- study_tips[nzchar(study_tips)]
unseen <- setdiff(study_tips, tr$tip.label)
if (length(unseen))
  warning("study tips absent from the tree: ", paste(unseen, collapse = ", "))

# A tip named <allele>_<sitecode>, e.g. blaCTX-M-36_S1W, is rendered as
# "blaCTX-M-36 (Sea 1, water)": a known allele recovered in this study, not a
# new haplotype. Anything else keeps its own name.
suffix_site <- function(lab) {
  m <- regmatches(lab, regexpr(paste0("_(", paste(SITEMAP$label, collapse = "|"),
                                      ")$"), lab))
  if (!length(m)) return(NA_character_)
  sub("^_", "", m)
}

study <- data.frame(label = study_tips, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(site_suffix = suffix_site(label)) %>%
  ungroup() %>%
  left_join(SITEMAP, by = "label") %>%
  left_join(SITEMAP %>% rename(site_suffix = label,
                               sfx_loc = location, sfx_mat = matrix),
            by = "site_suffix") %>%
  mutate(
    show = case_when(
      !is.na(location)    ~ sprintf("%s - %s, %s", label, location, matrix),
      !is.na(sfx_loc)     ~ sprintf("%s (%s, %s)",
                                    sub("_[^_]+$", "", label), sfx_loc, sfx_mat),
      TRUE                ~ label),
    class = case_when(
      is.na(location) & !is.na(sfx_loc) ~ "Known allele detected in this study",
      is.na(location)                   ~ "This study - additional sequence",
      matrix == "plastic"               ~ "This study - plastic",
      TRUE                              ~ "This study - water")) %>%
  select(label, show, class)
message("study tips: ", paste(study$label, collapse = ", "))

tips <- tr$tip.label
coll <- grepl("_and_", tips)
og   <- grepl(OGPAT, tips)

meta <- data.frame(label = tips, stringsAsFactors = FALSE) %>%
  left_join(study, by = "label") %>%
  mutate(
    class = case_when(!is.na(class)           ~ class,
                      og                      ~ "Group 9 outgroup",
                      label == "blaCTX-M-177" ~ "Closest reference (CTX-M-177)",
                      coll                    ~ "Collapsed reference haplotype",
                      TRUE                    ~ "Reference allele"),
    show  = case_when(!is.na(show) ~ show,
                      coll ~ sub("_and_(\\d+)_more", " (+\\1)", label),
                      TRUE ~ label),
    bold  = ifelse(is.na(match(label, study$label)), "plain", "bold"),
    expr  = gene_exprs(show, bold == "bold"))

# Support is shown at every node carrying a value; nodes with SH-aLRT >= 80 and
# UFBoot >= 95 are marked with a filled dot.
#
# IQ-TREE writes "SH-aLRT/UFBoot" into the internal-node labels, which ape keeps
# in $node.label. Node i of $node.label is tree node Ntip + i.
nlab   <- tr$node.label %||% character(0)
paired <- grepl("/", nlab)
a <- suppressWarnings(as.numeric(sub("/.*", "", nlab)))
u <- suppressWarnings(as.numeric(sub(".*/", "", nlab)))
a[!paired] <- NA; u[!paired] <- NA
supdf <- data.frame(
  node   = seq_along(nlab) + length(tr$tip.label),
  suplab = ifelse(!is.na(a) & !is.na(u), sprintf("%.0f/%.0f", a, u), NA_character_),
  ok     = !is.na(a) & !is.na(u) & a >= 80 & u >= 95)
message("support labels parsed: ", sum(!is.na(supdf$suplab)),
        " | passing >=80/>=95: ", sum(supdf$ok))

pal <- c("This study - water & plastic"  = "#C0392B",
         "This study - plastic"          = "#D35400",
         "This study - water"            = "#2980B9",
         "Closest reference (CTX-M-177)" = "#E67E22",
         "Collapsed reference haplotype" = "#7F8C8D",
         "Reference allele"              = "#2C3E50",
         "Group 9 outgroup"              = "#16A085",
         "Known allele detected in this study" = "#8E44AD",
         "This study - additional sequence"   = "#8E44AD")
shp <- c("This study - water & plastic" = 16, "This study - plastic" = 17,
         "This study - water" = 16, "Closest reference (CTX-M-177)" = 18,
         "Collapsed reference haplotype" = 15, "Reference allele" = 1,
         "Group 9 outgroup" = 8,
         "Known allele detected in this study" = 23,
         "This study - additional sequence"   = 15)

# show.legend = FALSE on the label layer: without it ggplot adds the text glyph
# to every legend key and each entry renders as a stray "a" beside its symbol.
p <- ggtree(tr, size = 0.35) %<+% meta %<+% supdf

# Tip labels are plotmath (italic bla, subscript allele), hence parse = TRUE.
# Support sits above the branch it belongs to, centred on it, so it no longer
# collides with the tip symbols.
p <- p +
  geom_tiplab(aes(label = expr, colour = class), parse = TRUE,
              size = 2.7, offset = 0.0012, align = TRUE,
              linetype = "dotted", linesize = 0.15, show.legend = FALSE) +
  geom_tippoint(aes(colour = class, shape = class), size = 1.8) +
  geom_text2(aes(subset = !is.na(suplab) & !isTip &
                   (is.null(STEM) | node != STEM$stem_child %||% -1), x = branch, label = suplab),
             size = 1.8, vjust = -0.5, colour = "grey30") +
  geom_nodepoint(aes(subset = ok), size = 1.1, colour = "black") +
  scale_colour_manual(values = pal, breaks = names(pal)) +
  scale_shape_manual(values  = shp, breaks = names(shp)) +
  geom_treescale(x = 0, y = -1.5, width = 0.005, fontsize = 2.6,
                 offset = 0.6, linesize = 0.4)

if (!is.null(STEM)) {
  d  <- p$data
  y  <- d$y[d$node == STEM$stem_child]
  x0 <- d$x[d$node == STEM$stem_parent]; x1 <- d$x[d$node == STEM$stem_child]
  xm <- x0 + (x1 - x0) * 0.3; dx <- (x1 - x0) * 0.04; dy <- 0.45
  brk <- data.frame(x = c(xm - 2 * dx, xm + dx) - dx, xend = c(xm, xm + 3 * dx) - dx,
                    y = y - dy, yend = y + dy)
  p <- p +
    annotate("rect", xmin = xm - 2 * dx, xmax = xm + dx, ymin = y - 0.3,
             ymax = y + 0.3, fill = "white", colour = NA) +
    geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, linewidth = 0.35) +
    annotate("text", x = x0 - (x1 - x0) * 0.3, y = y, size = 2.2, hjust = 1,
             label = sprintf("%.2f\n(not to scale)", STEM$true_len)) +
    annotate("text", x = x0 + (x1 - x0) * 0.68, y = y - 1.0, size = 1.8,
             colour = "grey30", label = STEM_SUP) +
    hexpand(0.08, direction = -1)
}

p <- p +
  theme_tree() +
  # legend below the panel: at 66 tips an inset legend covers the branches of
  # whichever clade it lands on, which is how the first version hid them.
  theme(legend.position   = "bottom",
        legend.title      = element_blank(),
        legend.text       = element_text(size = 8),
        legend.key.height = unit(0.38, "cm"),
        plot.margin       = margin(6, 6, 6, 6)) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE),
         shape  = guide_legend(nrow = 3, byrow = TRUE)) +
  hexpand(0.45) +                      # room for the longest label
  coord_cartesian(clip = "off")

pdf_out <- file.path(OUTDIR, "Figure_CTXM_sample_tree.pdf")
png_out <- file.path(OUTDIR, "Figure_CTXM_sample_tree.png")
ggsave(pdf_out, p, width = 7.5, height = 11, limitsize = FALSE)
ggsave(png_out, p, width = 7.5, height = 11, dpi = 600, limitsize = FALSE)
cat(pdf_out, " / ", png_out, " written\n", sep = "")
cat("Zero-length branches collapsed. SH-aLRT/UFBoot shown at every remaining node;\n",
    "filled dots mark >=80/>=95. Group 1 stem shortened for display.\n")
