# blaCTX-M group 1 amplicon phylogeny

Maximum-likelihood tree of the blaCTX-M group 1 amplicon consensus sequences
recovered in this study, placed among reference group 1 alleles and rooted on
group 9. This tree is separate from the 16S reference tree in `phylogeny/`.

```
study consensus sequences + curated group 1 references + group 9 outgroup
  -> MAFFT -> CDS positions 28-400 -> collapse identical references -> IQ-TREE
```

## Scripts

| script | what it does | output |
|---|---|---|
| `01_build_panel_and_tree.sh` | curates the reference panel, aligns, cuts to the window, collapses identical alleles, infers the tree | `output/ctxm/` |
| `02_plot_complete_tree.R` | complete tree with all tips and support values | `output/figures/Figure_CTXM_sample_tree.pdf/.png` |
| `03_plot_figure5c.R` | Figure 5C, the tree pruned to the tips in `figure5c_tips.txt` | `output/figures/Figure5C_CTXM.pdf/.png` |
| `plot_helpers.R` | functions shared by the two plotting scripts | |
| `run_all.sh` | runs all three in order | |

## Running

The tree is deposited in `output/ctxm/`, so the figures can be redrawn from the
repository alone. From the repository root:

```bash
Rscript phylogeny/ctxm/02_plot_complete_tree.R
Rscript phylogeny/ctxm/03_plot_figure5c.R
```

To rebuild the tree from the consensus and reference sequences:

```bash
export CTXM_DATA=/path/to/ctxm/analysis   # holds refs/ and consensus/
bash phylogeny/ctxm/run_all.sh
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
(`output/ctxm/drop.txt`):

| allele | reason |
|---|---|
| CTX-M-211 | 9 bp insertion that disrupts the CDS coordinates |
| CTX-M-123, CTX-M-132, CTX-M-153, CTX-M-234 | group 1/group 9 recombinants |

The remaining 135 references were collapsed to 58 tips, one per sequence that
is distinct over the window. `output/ctxm/haplotype_collapse.tsv` lists the
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

All values are from `output/ctxm/ctxm_group1.iqtree` and
`output/ctxm/ctxm_versions.txt`.

The six study sequences are the consensus sequences from five libraries (R1W,
R1P, R2W, S1W, ScW) and a second sequence from S1W that is identical to
blaCTX-M-36. The five consensus sequences are identical over the window. S1P
and ScP, with 39 and 17 reads, did not yield a consensus.

## Figures

Both figures are drawn from `output/ctxm/ctxm_group1.treefile`.

- Branches of zero length (38 of 61 internal branches) are collapsed to
  polytomies before drawing.
- Support is shown as SH-aLRT/UFBoot. Nodes with SH-aLRT >= 80 and UFBoot >= 95
  are considered well supported.
- Gene names are set as italic *bla* with the allele as subscript in the
  complete tree, and as bold italic labels in Figure 5C.

**Complete tree** (`02_plot_complete_tree.R`). All 64 tips, including the
group 9 outgroup. The branch joining group 1 to the outgroup (0.47
substitutions per site) is drawn shortened and marked with a break. Support is
shown at every node; filled dots mark nodes passing both thresholds.

**Figure 5C** (`03_plot_figure5c.R`). The rooted tree pruned to the tips in
`figure5c_tips.txt`; the outgroup is not drawn. Support at each node is that of
the same clade in the complete tree. The two branches leaving the root are
drawn at a fixed length and marked with a break. The filled star marks the
blaCTX-M-177-like sequence shared by R1W, R1P, R2W, S1W and ScW (drawn once,
as they are identical); the open star marks blaCTX-M-36, identical to the
second S1W sequence.

The panel used in the composite Figure 5 is the same script at a smaller size;
the command is in the top-level `README.md`.

## Files in output/ctxm

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
