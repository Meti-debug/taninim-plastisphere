#!/usr/bin/env bash
# =============================================================================
# CTX-M group 1 amplicon phylogeny.
#
# INPUT   $CTXM_DATA/refs/group1_v2.fasta        curated group-1 reference panel
#         $CTXM_DATA/refs/group9_v2.fasta        group-9 alleles, used as outgroup
#         $CTXM_DATA/consensus/*.consensus.fasta per-sample amplicon consensus
#
# OUTPUT  output/ctxm/group1_clean.fasta         panel after curation
#         output/ctxm/set.aln.fasta             MAFFT alignment
#         output/ctxm/tree_input.fasta           trimmed, collapsed, relabelled
#         output/ctxm/haplotype_collapse.tsv     which alleles each tip stands for
#         output/ctxm/ctxm_versions.txt          tool versions for the Methods
#         output/ctxm/ctxm_group1.treefile       the tree, Newick
#         output/ctxm/ctxm_group1.iqtree         model fit and alignment stats
#
# Run from the repository root:
#
#   export CTXM_DATA=/path/to/recovered_ctxm/analysis
#   bash phylogeny/ctxm/01_build_panel_and_tree.sh
#
# Requires mafft, iqtree2 and seqkit on PATH, plus biopython:
#   conda create -n phylo -c conda-forge -c bioconda iqtree mafft seqkit -y
#
# Unlike the 16S tree in phylogeny/, this alignment is NOT trimmed with trimAl.
# It is cut to a fixed coordinate window, because the amplicon covers only part
# of the CDS and the window has to be the same in every sample for the tips to
# be comparable. Report the window, not a trimming heuristic.
# =============================================================================
set -euo pipefail

D="${CTXM_DATA:?export CTXM_DATA to the directory holding refs/ and consensus/}"
WORK="${CTXM_WORK:-output/ctxm}"
LO="${CTXM_LO:-1}"          # first CDS position covered by the amplicon
HI="${CTXM_HI:-400}"        # last position with usable depth (>=100x, all samples)
THREADS="${CTXM_THREADS:-8}"

for f in "$D/refs/group1_v2.fasta" "$D/refs/group9_v2.fasta"; do
    [ -s "$f" ] || { echo "missing $f"; exit 1; }
done
ls "$D"/consensus/*.consensus.fasta >/dev/null 2>&1 || {
    echo "no consensus FASTA files in $D/consensus/"; exit 1; }

mkdir -p "$WORK"
cd "$WORK"

IQ=$(command -v iqtree2 || command -v iqtree) || {
    echo "iqtree not found on PATH"; exit 1; }

# --- 1. Curate the reference panel -------------------------------------------
# These five alleles are removed deliberately, and the reasons belong in the
# Methods and the supplementary, not only here:
#   CTX-M-211             885 bp, a 9 bp insertion that breaks CDS coordinates
#   CTX-M-123,-132        90-93% identity to the group-1 consensus; group-1/9
#   CTX-M-153,-234        recombinants carrying the same signature
printf 'blaCTX-M-%s\n' 123 132 153 234 211 > drop.txt
seqkit grep -v -n -f drop.txt "$D"/refs/group1_v2.fasta > group1_clean.fasta
echo "panel: $(grep -c '^>' group1_clean.fasta) group-1 alleles retained"

# --- 2. Outgroup: CTX-M group 9, a distinct lineage --------------------------
# Rooting on group 9 is an outgroup rooting. Rooting on a group-1 allele
# instead is an ingroup rooting and does not locate the root of group 1; see
# phylogeny/ctxm/README.md before changing this.
seqkit grep -nrp 'blaCTX-M-(14|27|65)$' "$D"/refs/group9_v2.fasta > outgroup.fasta
echo "outgroup: $(grep -c '^>' outgroup.fasta) group-9 alleles"

# --- 3. Record versions now, so the Methods cannot drift from the run --------
{
    echo "run date   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "mafft      $(mafft --version 2>&1 | head -1)"
    echo "iqtree     $($IQ --version 2>&1 | head -1)"
    echo "seqkit     $(seqkit version 2>&1 | head -1)"
    echo "window     CDS positions $LO-$HI"
    echo "model      MFP (ModelFinder; selection recorded in ctxm_group1.iqtree)"
    echo "support    1000 SH-aLRT replicates, 1000 ultrafast bootstrap, --bnni"
} > ctxm_versions.txt
cat ctxm_versions.txt

# --- 4. Assemble the study sequences -----------------------------------------
# One consensus per library, plus any additional study haplotypes named in
# CTXM_EXTRA. A single consensus per library collapses within-sample variation,
# so a second allele present in one library is invisible unless it is supplied
# here as its own record. Example:
#
#   export CTXM_EXTRA="$HOME/ctxm_rebuild/hap2_S1W.fasta"
#
# Study records are identified by membership in study_ids.txt rather than by an
# ID prefix, so extra files keep whatever headers they already carry.
cat "$D"/consensus/*.consensus.fasta > study.fasta
if [ -n "${CTXM_EXTRA:-}" ]; then
    for f in $CTXM_EXTRA; do
        [ -s "$f" ] || { echo "CTXM_EXTRA file not found or empty: $f"; exit 1; }
        cat "$f" >> study.fasta
        echo "extra study sequences from $f: $(grep -c '^>' "$f")"
    done
fi
grep '^>' study.fasta | sed 's/^>//;s/[[:space:]].*//' > study_ids.txt
echo "study sequences: $(wc -l < study_ids.txt)"

# --- 5. Align (--adjustdirection guards against reverse-complement input) ----
cat group1_clean.fasta outgroup.fasta study.fasta > set.fasta
mafft --auto --adjustdirection --thread "$THREADS" set.fasta > set.aln.fasta

# --- 6. Trim to the window, collapse indistinguishable alleles, relabel ------
python3 - "$LO" "$HI" <<'PYEOF'
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

LO, HI = int(sys.argv[1]) - 1, int(sys.argv[2])

# sample index -> transect site (W = water, P = plastic biofilm)
SITE = {"s01": "R1W", "s02": "R1P", "s03": "R2W", "s04": "S1W",
        "s05": "S1P", "s06": "ScW", "s07": "ScP"}

# mafft --adjustdirection prefixes reverse-complemented records with _R_, so
# strip that before matching an ID against the study list.
def base_id(i):
    return i[3:] if i.startswith("_R_") else i

study_ids = {l.strip() for l in open("study_ids.txt") if l.strip()}
recs = list(SeqIO.parse("set.aln.fasta", "fasta"))
qs   = [r for r in recs if base_id(r.id) in study_ids]
refs = [r for r in recs if base_id(r.id) not in study_ids]
missing = study_ids - {base_id(r.id) for r in qs}
if missing:
    sys.exit("study sequences lost during alignment: " + ", ".join(sorted(missing)))
if not qs:
    sys.exit("no study sequences found; check study_ids.txt")

# collapse references identical across the amplified window
groups = defaultdict(list)
for r in refs:
    groups[str(r.seq[LO:HI]).upper()].append(r.id)

out, table = [], []
for seq, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    members = sorted(members)
    if len(members) == 1:
        name = members[0]
    else:
        # name the collapsed tip after its most clinically recognised member
        pref = [m for m in members if m in
                ("blaCTX-M-15", "blaCTX-M-1", "blaCTX-M-3", "blaCTX-M-55")]
        head = pref[0] if pref else members[0]
        name = f"{head}_and_{len(members)-1}_more"
    out.append(SeqRecord(Seq(seq), id=name, description=""))
    table.append((name, len(members), ",".join(members)))

def tip_name(rid):
    """Per-library consensus -> site code; anything else keeps its own name."""
    rid = base_id(rid)
    parts = rid.split("_")
    if len(parts) > 1 and parts[1] in SITE:
        return SITE[parts[1]]
    return rid

study_tips = []
for r in qs:
    name = tip_name(r.id)
    study_tips.append(name)
    out.append(SeqRecord(Seq(str(r.seq[LO:HI]).upper()), id=name, description=""))

SeqIO.write(out, "tree_input.fasta", "fasta")
with open("haplotype_collapse.tsv", "w") as fh:
    fh.write("tip_label\tn_alleles\tmembers\n")
    for n, c, m in table:
        fh.write(f"{n}\t{c}\t{m}\n")
# the plotting scripts read this rather than carrying a hardcoded sample table,
# so an extra haplotype supplied through CTXM_EXTRA is coloured automatically
with open("study_tips.txt", "w") as fh:
    fh.write("\n".join(study_tips) + "\n")

print(f"{len(refs)} alleles -> {len(groups)} haplotypes; "
      f"{len(qs)} study sequences; {len(out)} tips")
print("study tips: " + ", ".join(study_tips))
print("collapse table -> haplotype_collapse.tsv  (cite this in the Methods)")
PYEOF

# --- 7. Maximum-likelihood tree ----------------------------------------------
# -m MFP  ModelFinder selects the model; do not hardcode the selected model
# -B      ultrafast bootstrap      -alrt  SH-aLRT
# --bnni  reduces UFBoot bias from model violation
# The collapse step in step 6 renames any tip that stands for more than one
# allele, so blaCTX-M-14 can become blaCTX-M-14_and_1_more. Read the outgroup
# names back from the tree input rather than from outgroup.fasta, or IQ-TREE
# exits with "Alignment does not have specified outgroup taxon". All group-9
# tips are passed together so the outgroup is rooted as a clade.
OGPAT=$(grep '^>' outgroup.fasta | sed 's/^>//;s/ .*//;s/^blaCTX-M-//' | paste -sd'|')
OG=$(grep '^>' tree_input.fasta | sed 's/^>//' \
     | grep -E "^blaCTX-M-(${OGPAT})(_and_[0-9]+_more)?$" | paste -sd,)
[ -n "$OG" ] || { echo "no group-9 outgroup tip found in tree_input.fasta"; exit 1; }
echo "outgroup tips: $OG"
"$IQ" -s tree_input.fasta \
      -m MFP -B 1000 -alrt 1000 --bnni \
      -T "$THREADS" -o "$OG" \
      --prefix ctxm_group1 -redo

[ -s ctxm_group1.treefile ] || { echo "no tree produced; see ctxm_group1.log"; exit 1; }

echo
grep -E 'Best-fit model|Number of parsimony informative sites' ctxm_group1.iqtree || true
echo
echo "tree      $WORK/ctxm_group1.treefile"
echo "report    $WORK/ctxm_group1.iqtree     (model, alignment length, support)"
echo "versions  $WORK/ctxm_versions.txt      (quote these in the Methods)"
echo "deposit   $WORK/haplotype_collapse.tsv"
echo
echo "next: Rscript phylogeny/ctxm/02_plot_tree_gene.R"
