# =============================================================================
# Shared helpers for 02_plot_complete_tree.R and 03_plot_figure5c.R.
#
#   source("analysis/phylogeny/ctxm/plot_helpers.R")
# =============================================================================

# Gene names as in the paper: italic "bla", subscript allele.
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
# what is needed to mark the break and label its true length.
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


# The tree drawn in Figure 5C and Figure S6, prepared identically for both:
#   1. root on the group 9 outgroup;
#   2. collapse zero-length branches in the complete tree;
#   3. prune to the tips listed in TIPS (the outgroup is not drawn);
#   4. take support for each drawn node from the same clade of the complete tree;
#   5. draw the two root branches at a fixed length ROOT_LEN.
panel_tree <- function(TREE, TIPS, ROOT_LEN, HAPNAME) {
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
  list(tree = tr, support = supdf, n = n, tip_label = tip_label)
}
