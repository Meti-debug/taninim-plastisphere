#!/usr/bin/env bash
# Build the reference 16S tree end to end.
#
#   bash phylogeny/run_all.sh
#
# Run from the repository root, with data/counts.csv in place and
# NCBI_EMAIL set.  Roughly 30-60 minutes: most of it is the NCBI retrieval and
# the IQ-TREE run.

set -euo pipefail
[ -s data/counts.csv ] || { echo "put counts.csv in data/ first"; exit 1; }
: "${NCBI_EMAIL:?export NCBI_EMAIL=you@example.org  (NCBI requires a contact address)}"
mkdir -p output

echo "== 1/2  retrieving reference sequences"
python3 phylogeny/01_fetch_reference_16S.py

echo "== 2/2  aligning and inferring the tree"
bash phylogeny/02_align_trim_tree.sh

echo
echo "tree      output/phylo_tree.treefile"
echo "rings     written by figures/fig4c_itol_annotations.R"
echo "deposit   output/phylo_tip_provenance.csv"
