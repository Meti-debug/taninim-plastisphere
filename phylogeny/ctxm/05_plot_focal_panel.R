#!/usr/bin/env Rscript
# =============================================================================
# CTX-M group 1 tree, reduced to a figure panel.
#
# The full tree has 70 tips and needs a page of its own. This keeps the study
# sequences, the group 9 outgroup, and the K reference alleles nearest any
# study sequence by patristic distance, which is reproducible and does not
# depend on eyeballing the full tree.
#
#   Rscript phylogeny/ctxm/05_plot_focal_panel.R
#   CTXM_KEEP=14 Rscript phylogeny/ctxm/05_plot_focal_panel.R   # wider context
#
# IMPORTANT: pruning a tree can misalign $node.label in some ape versions, so
# the script counts the support labels before and after and refuses to write a
# figure whose support values were lost. If it stops, report the counts.
# =============================================================================
suppressPackageStartupMessages({
  library(ggtree); library(ggplot2); library(ape); library(dplyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

TREE   <- Sys.getenv("CTXM_TREE",   "output/ctxm/ctxm_group1.treefile")
OUTDIR <- Sys.getenv("CTXM_OUTDIR", "figures")
KEEP_N <- as.integer(Sys.getenv("CTXM_KEEP", "10"))
LEGEND  <- Sys.getenv("CTXM_LEGEND", "bottom")   # "bottom" or "none"
SUPPORT <- Sys.getenv("CTXM_SUPPORT", "show")    # "show" or "none"
OFFSET  <- as.numeric(Sys.getenv("CTXM_OFFSET", "0.0006"))  # tip symbol to label
PANEL_OG <- Sys.getenv("CTXM_PANEL_OG", "yes")   # "no" drops the outgroup from
# the drawing. The root-to-outgroup branch is ~0.5 while the whole ingroup spans
# ~0.01, so at panel size the outgroup flattens the ingroup into a vertical
# line. Dropping it makes the panel legible; the caption then has to state the
# rooting, and the full tree shows it.
FIG_W  <- as.numeric(Sys.getenv("CTXM_W", "7.5"))
FIG_H  <- as.numeric(Sys.getenv("CTXM_H", "0"))  # 0 = scale with tip count
OGPAT  <- "^blaCTX-M-(14|27|65)(_and_[0-9]+_more)?$"

if (!file.exists(TREE)) stop("tree not found: ", TREE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

tr <- ape::read.tree(TREE)
n_before <- sum(grepl("/", tr$node.label %||% character(0)))

TIPFILE <- Sys.getenv("CTXM_STUDY_TIPS",
                      file.path(dirname(TREE), "study_tips.txt"))
if (!file.exists(TIPFILE)) stop("study_tips.txt not found at ", TIPFILE)
study_tips <- trimws(readLines(TIPFILE, warn = FALSE))
study_tips <- study_tips[nzchar(study_tips)]
og_tips    <- grep(OGPAT, tr$tip.label, value = TRUE)

# nearest reference alleles, by patristic distance to any study sequence
d <- ape::cophenetic.phylo(tr)
cand <- setdiff(tr$tip.label, c(study_tips, og_tips))
nearest <- sort(apply(d[cand, study_tips, drop = FALSE], 1, min))
keep <- unique(c(study_tips,
                 if (PANEL_OG == "no") character(0) else og_tips,
                 names(head(nearest, KEEP_N))))
message("keeping ", length(keep), " of ", length(tr$tip.label), " tips")

sub <- ape::keep.tip(tr, keep)

# Support on a pruned tree, done honestly.
#
# Pruning destroys clades. The one well-supported node in this tree is the
# group 1 / group 9 split at 100/100, and it spans all 68 group 1 tips; a panel
# keeping 16 of them does not draw that clade at all, so no support value on the
# panel can stand for it. Taking the support of each displayed clade's MRCA in
# the full tree is worse than nothing here, because that ancestor spans 56 tips
# where the drawn clade spans 16.
#
# A displayed node therefore inherits a label only when exactly one node of the
# full tree restricts to its tip set, meaning the drawn split really is that
# split. Where a path collapsed, the node is left blank. CTXM_SUPPORT=none omits
# the values entirely, which is the right choice at panel size: say in the
# caption that the tree is pruned for display and that support is shown on the
# full tree.
clade_sets <- function(phy) {
  pp <- ape::prop.part(phy)
  labs <- attr(pp, "labels")
  out <- lapply(pp, function(ix) sort(labs[ix]))
  if (length(out) != phy$Nnode)
    stop("prop.part returned ", length(out), " clades for ", phy$Nnode,
         " internal nodes; the node mapping cannot be trusted")
  if (length(out[[1]]) != length(phy$tip.label))
    stop("prop.part element 1 is not the root clade; node mapping unreliable")
  out
}

full_clades <- clade_sets(tr)
full_labels <- tr$node.label %||% rep("", length(full_clades))
keyed <- list()
for (i in seq_along(full_clades)) {
  k <- paste(sort(intersect(full_clades[[i]], keep)), collapse = "\u0001")
  if (nzchar(k)) keyed[[k]] <- c(keyed[[k]], full_labels[i])
}

sub_labels <- vapply(clade_sets(sub), function(cl) {
  v <- unique(keyed[[paste(cl, collapse = "\u0001")]] %||% character(0))
  v <- v[grepl("/", v)]
  if (length(v) == 1) v else NA_character_
}, character(1))

if (SUPPORT == "none") {
  message("support values omitted (CTXM_SUPPORT=none)")
  sub_labels[] <- NA_character_
} else {
  message("support matched unambiguously: ", sum(!is.na(sub_labels)), " of ",
          length(sub_labels), " internal nodes")
}

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

suffix_site <- function(lab) {
  m <- regmatches(lab, regexpr(paste0("_(", paste(SITEMAP$label, collapse = "|"),
                                      ")$"), lab))
  if (!length(m)) return(NA_character_)
  sub("^_", "", m)
}

tips <- sub$tip.label
meta <- data.frame(label = tips, stringsAsFactors = FALSE) %>%
  rowwise() %>% mutate(sfx = suffix_site(label)) %>% ungroup() %>%
  left_join(SITEMAP, by = "label") %>%
  left_join(SITEMAP %>% rename(sfx = label, sloc = location, smat = matrix),
            by = "sfx") %>%
  mutate(
    is_study = label %in% study_tips,
    show = case_when(
      !is.na(location) ~ sprintf("%s - %s, %s", label, location, matrix),
      !is.na(sloc)     ~ sprintf("%s (%s, %s)",
                                 sub("_[^_]+$", "", label), sloc, smat),
      grepl("_and_", label) ~ sub("_and_(\\d+)_more", " (+\\1)", label),
      TRUE ~ label),
    class = case_when(
      grepl(OGPAT, label)              ~ "Group 9 outgroup",
      !is.na(location) & matrix == "plastic" ~ "This study - plastic",
      !is.na(location)                 ~ "This study - water",
      !is.na(sloc)                     ~ "Known allele detected in this study",
      label == "blaCTX-M-177"          ~ "Closest reference (CTX-M-177)",
      grepl("_and_", label)            ~ "Collapsed reference haplotype",
      TRUE                             ~ "Reference allele"),
    bold = ifelse(is_study, "bold", "plain")) %>%
  select(label, show, class, bold)

paired <- grepl("/", sub_labels) & !is.na(sub_labels)
a <- suppressWarnings(as.numeric(sub("/.*", "", sub_labels)))
u <- suppressWarnings(as.numeric(sub(".*/", "", sub_labels)))
a[!paired] <- NA; u[!paired] <- NA
supdf <- data.frame(
  node   = seq_along(sub_labels) + length(sub$tip.label),
  suplab = ifelse(paired, sprintf("%.0f/%.0f", a, u), NA_character_),
  ok     = paired & !is.na(a) & !is.na(u) & a >= 80 & u >= 95)
message("nodes passing >=80/>=95: ", sum(supdf$ok))

pal <- c("This study - plastic"                = "#D35400",
         "This study - water"                  = "#2980B9",
         "Known allele detected in this study" = "#8E44AD",
         "Closest reference (CTX-M-177)"       = "#E67E22",
         "Collapsed reference haplotype"       = "#7F8C8D",
         "Reference allele"                    = "#2C3E50",
         "Group 9 outgroup"                    = "#16A085")
shp <- c("This study - plastic" = 17, "This study - water" = 16,
         "Known allele detected in this study" = 23,
         "Closest reference (CTX-M-177)" = 18,
         "Collapsed reference haplotype" = 15, "Reference allele" = 1,
         "Group 9 outgroup" = 8)

p <- ggtree(sub, size = 0.5) %<+% meta %<+% supdf +
  geom_tiplab(aes(label = show, colour = class, fontface = bold),
              size = 3.0, offset = OFFSET, align = TRUE,
              linetype = "dotted", linesize = 0.2, show.legend = FALSE) +
  geom_tippoint(aes(colour = class, shape = class), size = 2.6) +
  geom_text2(aes(subset = !is.na(suplab), label = suplab),
             size = 2.3, hjust = 1.1, vjust = -0.5, colour = "grey35") +
  geom_nodepoint(aes(subset = ok), size = 1.5, colour = "black") +
  scale_colour_manual(values = pal, breaks = names(pal)) +
  scale_shape_manual(values  = shp, breaks = names(shp)) +
  geom_treescale(width = 0.002, fontsize = 2.8, offset = 0.5) +
  theme_tree() +
  theme(legend.position   = if (LEGEND == "none") "none" else "bottom",
        legend.title      = element_blank(),
        legend.text       = element_text(size = 7.5),
        legend.key.height = unit(0.35, "cm"),
        plot.margin       = margin(6, 60, 6, 6)) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE),
         shape  = guide_legend(nrow = 3, byrow = TRUE)) +
  hexpand(1.0) +
  coord_cartesian(clip = "off")

# a figure panel needs a fixed shape; a standalone figure can grow with tips
H <- if (FIG_H > 0) FIG_H else
     max(4.2, 0.22 * length(keep) + if (LEGEND == "none") 0.4 else 1.6)
pdf_out <- file.path(OUTDIR, "Figure_CTXM_panel.pdf")
png_out <- file.path(OUTDIR, "Figure_CTXM_panel.png")
ggsave(pdf_out, p, width = FIG_W, height = H, limitsize = FALSE)
ggsave(png_out, p, width = FIG_W, height = H, dpi = 600, limitsize = FALSE)
cat(pdf_out, " / ", png_out, " written (", length(keep), " tips, ",
    FIG_W, " x ", round(H, 2), " in)\n", sep = "")
cat(if (SUPPORT == "none")
      "Support omitted; show it on the full tree instead.\n"
    else
      "Support shown only where the drawn split matches the full tree exactly.\n")
