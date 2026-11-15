# =============================================================================
# Shared helpers for 02_plot_complete_tree.R and 03_plot_figure5c.R.
#
#   source("phylogeny/ctxm/plot_helpers.R")
# =============================================================================

# Gene names in journal style: italic "bla", subscript allele.
#   "blaCTX-M-36"                 -> italic bla, subscript CTX-M-36
#   "blaCTX-M-15 group (n = 56)"  -> the same, followed by " group (n = 56)"
#   anything without a bla prefix -> plain text
# Returns a plotmath string, so the layer that draws it needs parse = TRUE.
gene_expr <- function(lab, bold = FALSE) {
  q <- function(s) paste0('"', gsub('"', '\\\\"', s), '"')
  txt <- if (bold) "bold" else "plain"
  ita <- if (bold) "bolditalic" else "italic"
  # split into text pieces and gene names; every blaCTX-M-<n> is typeset
  pos <- gregexpr("bla(CTX-M-[0-9]+)", lab)[[1]]
  if (pos[1] == -1) return(sprintf("%s(%s)", txt, q(lab)))
  len <- attr(pos, "match.length")
  parts <- character(0); last <- 1
  for (k in seq_along(pos)) {
    if (pos[k] > last)
      parts <- c(parts, sprintf("%s(%s)", txt, q(substr(lab, last, pos[k] - 1))))
    allele <- substr(lab, pos[k] + 3, pos[k] + len[k] - 1)
    parts <- c(parts, sprintf('%s("bla")[%s(%s)]', ita, txt, q(allele)))
    last <- pos[k] + len[k]
  }
  if (last <= nchar(lab))
    parts <- c(parts, sprintf("%s(%s)", txt, q(substr(lab, last, nchar(lab)))))
  if (length(parts) == 1) parts else sprintf("paste(%s)", paste(parts, collapse = ", "))
}
gene_exprs <- function(labs, bold = FALSE) {
  bold <- rep_len(bold, length(labs))
  vapply(seq_along(labs), function(i) gene_expr(labs[i], bold[i]), character(1))
}

# Branches IQ-TREE leaves at its minimum length (about 1e-6) carry no
# substitution. They are polytomies drawn as bifurcations, and the support value
# on such a branch does not describe a resolved clade. Collapse them.
collapse_zero <- function(tr, tol = 1e-5) {
  n_before <- tr$Nnode
  tr <- ape::di2multi(tr, tol = tol)
  message("zero-length branches collapsed: ", n_before - tr$Nnode,
          " (internal nodes ", n_before, " -> ", tr$Nnode, ")")
  tr
}

# The stem joining group 1 to the group 9 outgroup is about 0.47 substitutions
# per site, while all of group 1 spans about 0.01. Drawn to scale, group 1
# collapses into a vertical line. Shorten that one edge for display and return
# what is needed to mark the break. The true length must be stated on the figure
# and in the caption.
shorten_stem <- function(tr, og_tips, show_len) {
  n <- length(tr$tip.label)
  ingroup <- setdiff(tr$tip.label, og_tips)
  stem_child <- ape::getMRCA(tr, ingroup)
  e <- which(tr$edge[, 2] == stem_child)
  if (length(e) != 1) stop("could not find the group 1 stem edge")
  true_len <- tr$edge.length[e]
  tr$edge.length[e] <- show_len
  message(sprintf("group 1 stem drawn at %.4f; true length %.4f", show_len, true_len))
  list(tree = tr, stem_child = stem_child, stem_parent = tr$edge[e, 1],
       true_len = true_len)
}

# SH-aLRT/UFBoot of the split separating group 1 from the outgroup.
stem_support <- function(tr, og_tips) {
  ingroup <- setdiff(tr$tip.label, og_tips)
  node <- ape::getMRCA(tr, ingroup)
  tr$node.label[node - length(tr$tip.label)]
}
