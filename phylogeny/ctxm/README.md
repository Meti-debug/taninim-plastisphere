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
| publication figure, gene-first tip labels | `02_plot_tree_gene.R` |
| publication figure with location and abundance panels | `03_plot_tree_bars.R` |
| one tip per sample, support values at every node | `04_plot_sample_tree.R` |

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

**Over the positions called in every sample, the seven sample sequences are
identical.** The branching among the study tips reflects which positions were
masked as N in each consensus, not sequence divergence. Nothing in that part of
the topology can support a claim about which lineages moved along the transect,
and the caption has to say so. Reporting the sample tips as distinct haplotypes
would misrepresent them.

Two consensus sequences do differ once the window is widened, which is why
`04_plot_sample_tree.R` carries a haplotype A/B assignment. That assignment is
window-dependent, so it is only meaningful alongside the window actually used.

## Rooting: unresolved, and it has to be resolved before submission

The two routes in this folder root the tree differently, and the difference is
not cosmetic.

- `01_build_panel_and_tree.sh` passes `-o <group-9 allele>` to IQ-TREE. CTX-M
  group 9 is a separate lineage, so this is an **outgroup rooting** and it does
  locate the root of group 1.
- `02_plot_tree_gene.R` and `03_plot_tree_bars.R` re-root in R on
  `blaCTX-M-36___g1`. CTX-M-36 is a **group-1** allele, so this is an **ingroup
  rooting**. It orients the drawing but does not locate the root of group 1,
  and the basal alleles that then appear at the foot of the figure are an
  artefact of that choice rather than a finding.

The published figure came from the second route. Whichever is kept, the Methods
must describe that one and must not describe both. If the group-9 rooting is
adopted, note also that `04_plot_sample_tree.R` drops the outgroup tips *after*
rooting, so the drawn tree is rooted without displaying its own outgroup — a
reader cannot otherwise reconcile the figure with the Methods.

## Support values: display is not interpretation

IQ-TREE is run with `-alrt 1000 -B 1000`, giving paired SH-aLRT/UFBoot labels.
`04_plot_sample_tree.R` prints the pair at **every** node and marks only those
passing SH-aLRT >= 80 and UFBoot >= 95.

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
coordinate window rather than trimmed with trimAl. Three windows appear across
the exploratory runs: 1-400, 8-422 and 28-402. **The window stated in the
Methods must be the one that produced the deposited tree**, and the retained
alignment length is reported by IQ-TREE in the `.iqtree` file.

The 28-402 window excludes the primer footprints. If a window including them is
used, say so, because primer-derived positions are not independent evidence.

## To report in the Methods

Measured from the run of 2026-11-10, recorded in
`output/ctxm/ctxm_versions.txt` and `output/ctxm/ctxm_group1.iqtree`:

```
mafft              v7.526, --auto --adjustdirection
seqkit             v2.13.0
iqtree             2.0.7, -m MFP -B 1000 -alrt 1000 --bnni
window             CDS positions 1-400
alignment          69 sequences, 400 columns, 93 distinct patterns
sites              127 parsimony-informative, 35 singleton, 238 constant
samples            5 amplicon consensus sequences (R1W, R1P, R2W, S1W, ScW)
reference tips     64 after collapsing identical alleles over the window
rooting            group-9 outgroup (blaCTX-M-14_and_1_more, blaCTX-M-65)
support shown      all nodes; interpreted as supported at >=80 SH-aLRT / >=95 UFBoot
reference alleles  135 collapsed to 64 distinct tips over the window
mafft strategy     --auto resolved to L-INS-i
model selected     K2P+R2, by BIC (AICc agreed; AIC preferred TIM3+F+R2)
log-likelihood     -1653.504
tree length        0.774
```

Take every number from those two files, not from memory.

The model is **selected**, not assumed. Writing "inference was performed under
K2P+R2" describes an a priori choice that was not made and invites a reviewer
to ask how it was justified. Write that ModelFinder selected it under BIC.

## Two things the run reports that belong in the paper

IQ-TREE prints, during alignment parsing:

```
NOTE: ScW is identical to R1W but kept for subsequent analysis
```

The sea-control and upstream-river consensus sequences are identical across
the whole 400 bp window. A statement that the haplotype was detected at all
four sites is true as detection, but it is four detections of one sequence,
not four related sequences, and the text should not imply otherwise.

Five libraries yielded a consensus, not seven. S1P and ScP, at 39 and 17
reads, did not. Any figure showing seven sample tips does not match this
alignment, and the Methods must state five.

Note also that the community 16S tree in `phylogeny/` was built in a separate
environment, with `apt-get install mafft trimal iqtree` and no version pin. Its
versions are not the ones recorded here and need their own check.
