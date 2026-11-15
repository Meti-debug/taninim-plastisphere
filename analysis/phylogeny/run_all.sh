#!/usr/bin/env bash
# Build the reference 16S tree end-to-end.
#
#   bash analysis/phylogeny/run_all.sh
#
# Run from the repository root, with data/counts.csv in place and
# the IQ-TREE run.

set -euo pipefail
[ -s data/counts.csv ] || { echo "put counts.csv in data/ first"; exit 1; }
: "${NCBI_EMAIL:?export NCBI_EMAIL=metigizawmaranata@gmail.com  (NCBI requires a contact address)}"
mkdir -p output

echo "== 1/2  retrieving reference sequences"
python3 analysis/phylogeny/01_fetch_reference_16S.py

echo "== 2/2  aligning and inferring the tree"
bash analysis/phylogeny/02_align_trim_tree.sh

echo
echo "tree      output/phylo_tree.treefile"
echo "annotate  in iTOL (itol.embl.de)"
echo "deposit   output/phylo_tip_provenance.csv"
