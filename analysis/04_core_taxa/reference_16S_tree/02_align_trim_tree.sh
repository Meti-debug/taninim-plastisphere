#!/usr/bin/env bash
# Align, trim, and infer the maximum-likelihood tree.
#
# INPUT   analysis/04_core_taxa/output/phylo_sequences_raw.fasta      from 01_fetch_reference_16S.py
# OUTPUT  analysis/04_core_taxa/output/phylo_aligned.fasta            MAFFT alignment
#         analysis/04_core_taxa/output/phylo_trimmed.fasta            trimAl output
#         analysis/04_core_taxa/output/phylo_tree.treefile            the tree, Newick
#         analysis/04_core_taxa/output/phylo_tree.iqtree              model fit and alignment stats
#         analysis/04_core_taxa/output/phylo_tree.log                 full IQ-TREE log
#         analysis/04_core_taxa/output/phylo_versions.txt             tool versions for the Methods
#
# Run from the repository root:  bash analysis/04_core_taxa/reference_16S_tree/02_align_trim_tree.sh
#
# Requires mafft, trimal and iqtree (or iqtree2) on PATH:
#   apt-get install -y mafft trimal iqtree
#
# Runtime for ~500 sequences: MAFFT 5-10 min, IQ-TREE 30-60 min on 4 threads.

set -euo pipefail

RAW=analysis/04_core_taxa/output/phylo_sequences_raw.fasta
ALN=analysis/04_core_taxa/output/phylo_aligned.fasta
TRIM=analysis/04_core_taxa/output/phylo_trimmed.fasta
PREFIX=analysis/04_core_taxa/output/phylo_tree
THREADS=AUTO               # IQ-TREE threads
MTHREADS=8                 # MAFFT threads
MODEL=GTR+G                  # model used for the reference tree
BOOTSTRAP=1000               # ultrafast bootstrap

[ -s "$RAW" ] || { echo "missing $RAW - run 01_fetch_reference_16S.py first"; exit 1; }

IQ=$(command -v iqtree2 || command -v iqtree) || {
    echo "iqtree not found on PATH"; exit 1; }

# ---- record versions now, so the Methods cannot drift from the run ---------
{
    echo "run date   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "mafft      $(mafft --version 2>&1 | head -1)"
    echo "iqtree     $($IQ --version 2>&1 | head -1)"
    echo "model      $MODEL"
    echo "bootstrap  $BOOTSTRAP ultrafast replicates"
} > analysis/04_core_taxa/output/phylo_versions.txt
cat analysis/04_core_taxa/output/phylo_versions.txt

# ---- align ------------------------------------------------------------------
echo "== MAFFT"
mafft --auto --reorder --thread "$MTHREADS" "$RAW" > "$ALN"
grep -c '^>' "$ALN" | xargs echo "  sequences aligned:"

# ---- trim -------------------------------------------------------------------
# -automated1 chooses between gappyout and strict from the alignment's own
# characteristics.  trim_alignment.py uses the trimal binary when it is on PATH
# and pytrimal otherwise (trimal is not packaged on Ubuntu 24.04), and FAILS
# rather than passing the untrimmed alignment forward: the original notebook
# fell through silently, which would leave the Methods describing a trimming
# step that never ran.
echo "== trimAl (automated1)"
python3 analysis/04_core_taxa/reference_16S_tree/trim_alignment.py "$ALN" "$TRIM"
[ -s "$TRIM" ] || { echo "trimming produced no output"; exit 1; }

# ---- infer ------------------------------------------------------------------
echo "== IQ-TREE"
"$IQ" -s "$TRIM" -m "$MODEL" -B "$BOOTSTRAP" -T "$THREADS" --prefix "$PREFIX" -redo

[ -s "${PREFIX}.treefile" ] || { echo "no tree produced - see ${PREFIX}.log"; exit 1; }
echo "tree: ${PREFIX}.treefile"
echo "next: open analysis/04_core_taxa/output/phylo_tree.treefile in iTOL (itol.embl.de)"
