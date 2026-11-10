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
CTXM_TREE="$WORK/ctxm_group1.treefile" Rscript phylogeny/ctxm/04_plot_sample_tree.R

# 02 and 03 expect the label format of the earlier tree behind the published
# figure, not the output of 01; see README.md. Run them explicitly with
# CTXM_TREE pointing at that tree.

echo
echo "tree      $WORK/ctxm_group1.treefile"
echo "report    $WORK/ctxm_group1.iqtree"
echo "versions  $WORK/ctxm_versions.txt"
echo "deposit   $WORK/haplotype_collapse.tsv"
echo "figures   figures/Figure_CTXM_sample_tree.pdf"
