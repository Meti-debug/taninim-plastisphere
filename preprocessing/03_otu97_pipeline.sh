#!/usr/bin/env bash
# 97% OTU table from full-length 16S Nanopore reads, all libraries clustered
# together in one pass.
#
#   IN_DIR=path/to/fastq  bash otu97_pipeline.sh
#
# INPUT   IN_DIR/<LIBRARY>.fastq or .fastq.gz   one file per library, named by
#         library (R1P1, R1W1, ...). Library names must be letters and digits
#         only: VSEARCH reads the library from the part of each read label
#         before the first dot.
# OUTPUT  OUT_DIR/otus.fasta           OTU representative sequences
#         OUT_DIR/otu_table.tsv        OTU x library read counts
#         OUT_DIR/read_tracking.tsv    reads per library after every step
#         OUT_DIR/run_settings.txt     VSEARCH version and every parameter used
#
# WHY IT IS BUILT THIS WAY
#   Nanopore reads carry ~1-3% error per base. Two reads of the same 16S gene
#   can therefore differ by close to 3%, the OTU threshold itself, and
#   clustering raw reads directly splits one organism into many OTUs. The
#   pipeline avoids that in three ways:
#     - OTUs are represented by CONSENSUS sequences (step 3), in which the
#       random errors of individual reads cancel out, not by single reads.
#     - All libraries are clustered TOGETHER, so the same organism gets the
#       same OTU in every library.
#     - Clusters are merged in order of ABUNDANCE (step 4), so the best-
#       supported sequences seed the final OTUs.
#   The OTUs are built from a fixed-size subsample of every library, which
#   keeps clustering tractable and gives each library equal weight. ALL
#   filtered reads are then counted against those OTUs (step 6).
#
# STEPS
#   1 filter     length and expected-error filter, reads labelled by library
#   2 subsample  BUILD_N reads per library, pooled
#   3 pass 1     cluster pooled reads at 97%, write one consensus per cluster,
#                drop clusters of a single read
#   4 pass 2     re-cluster the consensus sequences at 97%, largest first
#   5 chimeras   de novo chimera removal (UCHIME), relabel OTU_1, OTU_2, ...
#   6 map        every filtered read against the OTUs at 97% -> OTU table
#
# Each step is skipped if its output already exists, so an interrupted run
# can be restarted with the same command. Delete OUT_DIR to start over.

set -euo pipefail

# ---------------------------------------------------------------- settings
IN_DIR=${IN_DIR:?set IN_DIR to the folder of per-library FASTQ files}
OUT_DIR=${OUT_DIR:-otu97}
THREADS=${THREADS:-24}
VSEARCH=${VSEARCH:-vsearch}

MINLEN=${MINLEN:-1300}          # full-length 16S amplicon is ~1.45-1.55 kb
MAXLEN=${MAXLEN:-1700}
MAXEE_RATE=${MAXEE_RATE:-0.02}  # max expected errors per base (0.02 = 2%)
BUILD_N=${BUILD_N:-10000}       # reads per library used to build the OTUs
MIN_CLUSTER=${MIN_CLUSTER:-2}   # pass-1 clusters smaller than this are dropped
ID=${ID:-0.97}
SEED=${SEED:-11}

mkdir -p "$OUT_DIR"/{filtered,sub,work}
F="$OUT_DIR/filtered"; SUB="$OUT_DIR/sub"; W="$OUT_DIR/work"

shopt -s nullglob
FILES=("$IN_DIR"/*.fastq "$IN_DIR"/*.fastq.gz)
[[ ${#FILES[@]} -gt 0 ]] || { echo "no .fastq or .fastq.gz files in $IN_DIR" >&2; exit 1; }

lib_name() { local b; b=$(basename "$1"); b=${b%.gz}; echo "${b%.fastq}"; }
nseq() { grep -c '^>' "$1" || true; }

for fq in "${FILES[@]}"; do
    lib=$(lib_name "$fq")
    [[ $lib =~ ^[A-Za-z0-9]+$ ]] || {
        echo "library name '$lib' must be letters and digits only" >&2; exit 1; }
done

{
    echo "date        $(date -Iseconds)"
    echo "vsearch     $("$VSEARCH" --version 2>&1 | head -1)"
    echo "IN_DIR      $IN_DIR"
    echo "libraries   ${#FILES[@]}"
    for v in MINLEN MAXLEN MAXEE_RATE BUILD_N MIN_CLUSTER ID SEED THREADS; do
        printf '%-11s %s\n' "$v" "${!v}"
    done
} > "$OUT_DIR/run_settings.txt"
cat "$OUT_DIR/run_settings.txt"

# ---------------------------------------------------------------- 1 filter
# --fastq_qmax 93: Dorado reports qualities above 41, which VSEARCH otherwise
# rejects as invalid. --relabel R1P1. gives reads R1P1.1, R1P1.2, ...
echo "== 1 filter"
for fq in "${FILES[@]}"; do
    lib=$(lib_name "$fq")
    [[ -s "$F/$lib.fa" ]] && continue
    "$VSEARCH" --fastq_filter "$fq" \
        --fastq_qmax 93 \
        --fastq_minlen "$MINLEN" --fastq_maxlen "$MAXLEN" \
        --fastq_maxee_rate "$MAXEE_RATE" \
        --relabel "$lib." \
        --fastaout "$F/$lib.fa.tmp" --quiet
    mv "$F/$lib.fa.tmp" "$F/$lib.fa"
    echo "   $lib  $(nseq "$F/$lib.fa") reads kept"
done

# ---------------------------------------------------------------- 2 subsample
echo "== 2 subsample $BUILD_N reads per library"
for fa in "$F"/*.fa; do
    lib=$(basename "$fa" .fa)
    [[ -s "$SUB/$lib.fa" ]] && continue
    if (( $(nseq "$fa") <= BUILD_N )); then
        cp "$fa" "$SUB/$lib.fa"          # fewer reads than BUILD_N: use all
    else
        "$VSEARCH" --fastx_subsample "$fa" --sample_size "$BUILD_N" \
            --randseed "$SEED" --fastaout "$SUB/$lib.fa" --quiet
    fi
done
[[ -s "$W/pool.fa" ]] || cat "$SUB"/*.fa > "$W/pool.fa"
echo "   pooled $(nseq "$W/pool.fa") reads"

# ---------------------------------------------------------------- 3 pass 1
# --consout writes the majority-rule consensus of each cluster's aligned
# members; --sizeout records how many reads it represents.
echo "== 3 cluster pass 1 (consensus)"
if [[ ! -s "$W/pass1_kept.fa" ]]; then
    "$VSEARCH" --cluster_fast "$W/pool.fa" --id "$ID" --strand both \
        --threads "$THREADS" --sizeout \
        --consout "$W/pass1_cons.fa"
    "$VSEARCH" --sortbysize "$W/pass1_cons.fa" --sizein --sizeout \
        --minsize "$MIN_CLUSTER" --output "$W/pass1_kept.fa" --quiet
fi
echo "   $(nseq "$W/pass1_cons.fa") clusters, $(nseq "$W/pass1_kept.fa") with >= $MIN_CLUSTER reads"

# ---------------------------------------------------------------- 4 pass 2
echo "== 4 cluster pass 2 (consensus sequences, largest first)"
if [[ ! -s "$W/pass2.fa" ]]; then
    "$VSEARCH" --cluster_size "$W/pass1_kept.fa" --id "$ID" --strand both \
        --threads "$THREADS" --sizein --sizeout \
        --centroids "$W/pass2.fa"
fi
echo "   $(nseq "$W/pass2.fa") OTU candidates"

# ---------------------------------------------------------------- 5 chimeras
echo "== 5 chimera removal"
if [[ ! -s "$OUT_DIR/otus.fasta" ]]; then
    "$VSEARCH" --uchime_denovo "$W/pass2.fa" --sizein --sizeout \
        --nonchimeras "$W/nonchimeras.fa" --chimeras "$W/chimeras.fa" --quiet
    "$VSEARCH" --sortbysize "$W/nonchimeras.fa" --sizein \
        --relabel OTU_ --xsize --output "$OUT_DIR/otus.fasta" --quiet
fi
echo "   $(nseq "$W/chimeras.fa") chimeras removed, $(nseq "$OUT_DIR/otus.fasta") OTUs"

# ---------------------------------------------------------------- 6 map
# Every filtered read, not only the subsample, is assigned to its best OTU
# at >= 97% identity. Reads matching no OTU are counted as unassigned.
echo "== 6 map all filtered reads to OTUs"
if [[ ! -s "$OUT_DIR/otu_table.tsv" ]]; then
    cat "$F"/*.fa > "$W/all_filtered.fa"
    "$VSEARCH" --usearch_global "$W/all_filtered.fa" --db "$OUT_DIR/otus.fasta" \
        --id "$ID" --strand both --threads "$THREADS" \
        --otutabout "$OUT_DIR/otu_table.tsv.tmp"
    mv "$OUT_DIR/otu_table.tsv.tmp" "$OUT_DIR/otu_table.tsv"
    rm "$W/all_filtered.fa"
fi

# ---------------------------------------------------------------- tracking
# input reads are counted from the original FASTQ (4 lines per read)
echo "== read tracking"
awk -F'\t' 'NR==1 {for (i=2;i<=NF;i++) name[i]=$i; next}
            {for (i=2;i<=NF;i++) s[name[i]]+=$i}
            END {for (k in s) print k"\t"s[k]}' "$OUT_DIR/otu_table.tsv" \
    | sort > "$W/mapped_counts.tsv"
{
    printf 'library\tinput\tfiltered\tsubsampled\tmapped\tmapped_pct_of_filtered\n'
    for fq in "${FILES[@]}"; do
        lib=$(lib_name "$fq")
        if [[ $fq == *.gz ]]; then n_in=$(( $(zcat "$fq" | wc -l) / 4 ))
        else n_in=$(( $(wc -l < "$fq") / 4 )); fi
        n_f=$(nseq "$F/$lib.fa"); n_s=$(nseq "$SUB/$lib.fa")
        n_m=$(awk -v l="$lib" '$1==l {print $2}' "$W/mapped_counts.tsv"); n_m=${n_m:-0}
        pct=$(awk -v m="$n_m" -v f="$n_f" 'BEGIN {printf "%.1f", (f ? 100*m/f : 0)}')
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$lib" "$n_in" "$n_f" "$n_s" "$n_m" "$pct"
    done
} > "$OUT_DIR/read_tracking.tsv"
if command -v column >/dev/null; then column -t "$OUT_DIR/read_tracking.tsv"
else cat "$OUT_DIR/read_tracking.tsv"; fi
echo "done: $OUT_DIR/otu_table.tsv"
