# CTX-M group 1 amplicon phylogeny

    consensus/*.consensus.fasta + curated group-1 panel + group-9 outgroup
      -> MAFFT -> fixed coordinate window -> collapse identical refs -> IQ-TREE

This is a different tree from the one in `phylogeny/`. That one is built from
reference 16S sequences for the taxa detected in this study. This one is built
from the **amplicon consensus sequences of the samples themselves**, together
with reference alleles for context.

Run from the repository root:

```bash
export CTXM_DATA=/path/to/recovered_ctxm/analysis
bash phylogeny/ctxm/run_all.sh
```

Or one step at a time:

| step | script |
|---|---|
| curate the panel, align, trim to the window, infer the tree | `01_build_panel_and_tree.sh` |
| complete tree, one tip per sample, support at every resolved node | `04_plot_sample_tree.R` |
| Figure 5C: pruned panel, identical study tips drawn once | `06_plot_figure5c.R` |

Both plotting scripts source `plot_helpers.R` (gene-name typesetting,
zero-length branch collapse, stem shortening). Outputs go to `figures/`:
`Figure_CTXM_sample_tree.*` and `Figure5C_CTXM.*`.

Dependencies:

```bash
conda create -n phylo  -c conda-forge -c bioconda iqtree mafft seqkit -y
conda create -n rphylo -c conda-forge r-base r-ggplot2 r-dplyr r-ape \
                                     bioconductor-ggtree bioconductor-treeio -y
pip install biopython
```

`treeio` and `ggtree` are Bioconductor packages. `install.packages()` will not
find them and reports "not available for this version of R", which is
misleading — the problem is the repository, not the R version. Use
`BiocManager::install()` or the conda env above.

## What this tree is, and what it is not

The tips for this study are **per-sample consensus sequences of a single
blaCTX-M group 1 amplicon**, roughly 400 bp of the CDS. They are not whole
genes, not assemblies, and not separate alleles.

**Over the 28-400 window, the five sample consensus sequences are identical.**
Any branching among the study tips in the complete tree is IQ-TREE placing
identical sequences, not sequence divergence, and it cannot support a claim
about which lineages moved along the transect. Figure 5C therefore draws them as
one tip. Reporting the sample tips as distinct haplotypes would misrepresent
them.

## Rooting

`01_build_panel_and_tree.sh` passes `-o <group-9 allele>` to IQ-TREE. CTX-M
group 9 is a separate lineage, so this is an outgroup rooting and it locates the
root of group 1. Both figures draw the outgroup. Its stem is about 0.47
substitutions per site against about 0.01 for all of group 1, so that single
edge is shortened for display and marked with a break and its true length. In
the pruned Figure 5C the stem runs to the common ancestor of the retained group
1 tips and reads 0.51.

## Support values: display is not interpretation

IQ-TREE is run with `-alrt 1000 -B 1000`, giving paired SH-aLRT/UFBoot labels.
38 of the 61 internal branches are at IQ-TREE's minimum length and carry no
substitution. They are polytomies drawn as bifurcations, so both figures collapse
them first. `04_plot_sample_tree.R` then prints the pair at every remaining node
(23) and marks those passing SH-aLRT >= 80 and UFBoot >= 95. Figure 5C shows only
the group 1 versus group 9 split (100/100).

Low values therefore appear in the figure on purpose. That is not a
contradiction of the threshold as long as the Methods separates the two
statements: values are shown at all nodes, and nodes are regarded as well
supported at >= 80 / >= 95. Writing the threshold as a display rule makes the
figure look wrong.

## Reference curation

`01` drops five group-1 alleles before aligning, and the reasons belong in the
supplementary rather than only in a shell comment:

| allele | reason |
|---|---|
| CTX-M-211 | 885 bp; a 9 bp insertion breaks the CDS coordinates |
| CTX-M-123, -132 | 90-93% identity to the group-1 consensus; group-1/9 recombinants |
| CTX-M-153, -234 | same recombinant signature |

`01` also collapses reference alleles that are **identical over the amplified
window** into one tip, and writes `haplotype_collapse.tsv` listing which
alleles each tip stands for. **Deposit that table.** Without it a tip labelled
`blaCTX-M-15_and_11_more` cannot be interpreted, and the allele count in the
Methods cannot be checked against the panel.

## The window

The amplicon covers part of the CDS, so the alignment is cut to a fixed
coordinate window rather than trimmed with trimAl. The deposited tree uses CDS
positions 28-400, which excludes the primer footprints; this is the default in
`01_build_panel_and_tree.sh` (`CTXM_LO`, `CTXM_HI`). Exploratory runs used other
windows; if the window is changed, the Methods and every number below change
with it.

## To report in the Methods

Measured from the run of 2026-11-10, recorded in
`output/ctxm/ctxm_versions.txt` and `output/ctxm/ctxm_group1.iqtree`:

```
mafft              v7.526, --auto --adjustdirection
seqkit             v2.13.0
iqtree             2.0.7, -m MFP -B 1000 -alrt 1000 --bnni
window             CDS positions 28-400 (primer footprints excluded)
alignment          64 sequences, 373 columns, 82 distinct patterns
sites              116 parsimony-informative, 225 constant
study sequences    6: five consensus (R1W, R1P, R2W, S1W, ScW) and the
                   blaCTX-M-36-identical sequence from S1W
reference tips     58 after collapsing identical alleles over the window
reference alleles  135 collapsed to 58 distinct tips
rooting            group-9 outgroup (blaCTX-M-14_and_1_more, blaCTX-M-65)
model selected     K2P+R2, by BIC
log-likelihood     -1510.452
tree length        0.745
support shown      every non-zero-length node; supported at >=80 SH-aLRT / >=95 UFBoot
```

Take every number from those two files, not from memory.

The model is **selected**, not assumed. Write that ModelFinder selected it under
BIC.

## Two things the run reports that belong in the paper

IQ-TREE prints, during alignment parsing:

```
NOTE: R1P is identical to R1W but kept for subsequent analysis
NOTE: R2W (identical to R1W) is ignored but added at the end
NOTE: S1W (identical to R1W) is ignored but added at the end
NOTE: ScW (identical to R1W) is ignored but added at the end
NOTE: blaCTX-M-36_S1W is identical to blaCTX-M-36 but kept for subsequent analysis
```

The five consensus sequences are identical across the 373-column window, which
is why Figure 5C draws them as one tip. A statement that the haplotype was detected at all
four sites is true as detection, but it is four detections of one sequence,
not four related sequences, and the text should not imply otherwise.

Five libraries yielded a consensus, not seven. S1P and ScP, at 39 and 17
reads, did not.

Note also that the community 16S tree in `phylogeny/` was built in a separate
environment, with `apt-get install mafft trimal iqtree` and no version pin. Its
versions are not the ones recorded here and need their own check.
