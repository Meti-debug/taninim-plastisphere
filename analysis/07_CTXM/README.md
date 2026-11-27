# blaCTX-M group 1 amplicon phylogeny

Maximum-likelihood tree of the blaCTX-M group 1 amplicon consensus sequences
recovered in this study, placed among reference group 1 alleles and rooted on
group 9. This tree is separate from the 16S reference tree in `analysis/04_core_taxa/reference_16S_tree/`.

```
study consensus sequences + curated group 1 references + group 9 outgroup
  -> MAFFT -> CDS positions 28-400 -> collapse identical references -> IQ-TREE
```

Generates:
- Figure 5C
- Figure S6 A and B
- the complete CTX-M tree (supplementary)

Input:
- analysis/07_CTXM/tree/ (the deposited tree and alignment; enough to redraw the figures)
- $CTXM_DATA/refs and $CTXM_DATA/consensus (only to rebuild the tree)

Output:
- analysis/07_CTXM/output/Figure5C_CTXM.pdf/.png
- analysis/07_CTXM/output/FigureS6_AB.pdf/.png
- analysis/07_CTXM/output/Figure_CTXM_sample_tree.pdf/.png

`tree/` holds the deposited result and is tracked in git; `output/` holds
redrawn figures and is not.

## Scripts

| script | what it does | output |
|---|---|---|
| `01_build_panel_and_tree.sh` | curates the reference panel, aligns, cuts to the window, collapses identical alleles, infers the tree | `tree/` |
| `02_plot_complete_tree.R` | complete tree with all tips and support values | `output/Figure_CTXM_sample_tree.pdf/.png` |
| `03_plot_figure5c.R` | Figure 5C, the tree pruned to the tips in `figure5c_tips.txt` | `output/Figure5C_CTXM.pdf/.png` |
| `04_plot_figure_s6.R` | Figure S6 panels A and B: the same tree, with the variable sites read from `tree_input.fasta` | `output/FigureS6_AB.pdf/.png` |
| `plot_helpers.R` | functions shared by the plotting scripts, including the tree preparation used by both 03 and 04 | |

## Running

The tree is deposited in `analysis/07_CTXM/tree/`, so the figures can be redrawn from the
repository alone. From the repository root:

```bash
Rscript analysis/07_CTXM/02_plot_complete_tree.R
Rscript analysis/07_CTXM/03_plot_figure5c.R
```

To rebuild the tree from the consensus and reference sequences:

```bash
export CTXM_DATA=/path/to/ctxm/analysis   # holds refs/ and consensus/
bash analysis/07_CTXM/01_build_panel_and_tree.sh
```

## Dependencies

```bash
conda create -n phylo  -c conda-forge -c bioconda iqtree=2.0.7 mafft=7.526 seqkit -y
conda create -n rphylo -c conda-forge -c bioconda r-base r-ggplot2 r-dplyr r-ape \
                                     bioconductor-ggtree bioconductor-treeio -y
pip install biopython
```

`ggtree` and `treeio` are Bioconductor packages; outside conda, install them
with `BiocManager::install()`.

## Reference panel

Group 1 reference alleles were taken from the NCBI Bacterial Antimicrobial
Resistance Reference Gene Database. Five were excluded before alignment
(`analysis/07_CTXM/tree/drop.txt`):

| allele | reason |
|---|---|
| CTX-M-211 | 9 bp insertion that disrupts the CDS coordinates |
| CTX-M-123, CTX-M-132, CTX-M-153, CTX-M-234 | group 1/group 9 recombinants |

The remaining 135 references were collapsed to 58 tips, one per sequence that
is distinct over the window. `analysis/07_CTXM/tree/haplotype_collapse.tsv` lists the
alleles each tip represents. Tips are labelled `blaCTX-M-15_and_55_more`; the
figures show these as `blaCTX-M-15 group (n = 56)`, or `blaCTX-M-12 / -96` for
pairs.

## Alignment and tree

| | |
|---|---|
| alignment | MAFFT v7.526, `--auto --adjustdirection` |
| window | CDS positions 28-400, excluding the primer footprints |
| sequences | 64: 58 reference tips, 6 study sequences |
| columns | 373; 116 parsimony-informative, 225 constant |
| inference | IQ-TREE 2.0.7, `-m MFP -B 1000 -alrt 1000 --bnni` |
| model | K2P+R2, selected by ModelFinder under BIC |
| rooting | group 9 outgroup (CTX-M-14/-27, CTX-M-65) |
| log-likelihood | -1510.452 |

All values are from `analysis/07_CTXM/tree/ctxm_group1.iqtree` and
`analysis/07_CTXM/tree/ctxm_versions.txt`.

## Files in tree/

| file | content |
|---|---|
| `ctxm_group1.treefile` | maximum-likelihood tree with SH-aLRT/UFBoot labels |
| `ctxm_group1.contree` | UFBoot consensus tree |
| `ctxm_group1.iqtree`, `ctxm_group1.log` | IQ-TREE report and log |
| `ctxm_versions.txt` | tool versions, window and options |
| `tree_input.fasta` | the 373-column alignment used for the tree |
| `set.aln.fasta` | full-length alignment before cutting to the window |
| `group1_clean.fasta`, `outgroup.fasta`, `study.fasta` | reference, outgroup and study sequences |
| `haplotype_collapse.tsv` | alleles represented by each reference tip |
| `drop.txt` | alleles excluded before alignment |
| `study_ids.txt`, `study_tips.txt` | study sequence identifiers and tip labels |
