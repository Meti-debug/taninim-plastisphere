#!/usr/bin/env bash
# Build the CTX-M group 1 amplicon tree and the figures, end to end.
#
#   export CTXM_DATA=/path/to/recovered_ctxm/analysis
#   bash phylogeny/ctxm/run_all.sh
#
# Run from the repository root. A few minutes: the alignment is small and the
# tree has tens of tips, not hundreds.

set -euo pipefail
: "${CTXM_DATA:?export CTXM_DATA to the directory holding refs/ and consensus/}"

WORK="${CTXM_WORK:-output/ctxm}"
mkdir -p "$WORK" figures

echo "== 1/2  curating the panel, aligning, inferring the tree"
bash phylogeny/ctxm/01_build_panel_and_tree.sh

echo "== 2/2  figures"
CTXM_TREE="$WORK/ctxm_group1.treefile" Rscript phylogeny/ctxm/02_plot_complete_tree.R
CTXM_TREE="$WORK/ctxm_group1.treefile" Rscript phylogeny/ctxm/03_plot_figure5c.R

echo
echo "tree      $WORK/ctxm_group1.treefile"
echo "report    $WORK/ctxm_group1.iqtree"
echo "versions  $WORK/ctxm_versions.txt"
echo "deposit   $WORK/haplotype_collapse.tsv"
echo "figures   figures/Figure_CTXM_sample_tree.pdf  (complete tree)"
echo "          figures/Figure5C_CTXM.pdf            (Figure 5C)"
