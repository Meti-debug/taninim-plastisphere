# Phylogeny of the detected taxa

    data/counts.csv -> species list -> NCBI 16S -> MAFFT -> trimAl -> IQ-TREE

Everything at once, from the repository root:

```bash
export NCBI_EMAIL=you@example.org
bash phylogeny/run_all.sh
```

Or one step at a time:

| step | script |
|---|---|
| build the species list and fetch one 16S sequence per species | `01_fetch_reference_16S.py` |
| align, trim, infer the maximum-likelihood tree | `02_align_trim_tree.sh` |

The iTOL annotation files for the published figure are written by the R figure
code, not from this folder.

`trim_alignment.py` is a helper: it runs the `trimal` binary when it is on PATH
and pytrimal otherwise, because trimal is not packaged on Ubuntu 24.04.

Dependencies:

```bash
sudo apt-get install -y mafft iqtree
pip install biopython pandas numpy pytrimal
```

## What this tree is, and what it is not

The tips are **reference 16S sequences downloaded from NCBI for each species
name detected in this study**. They are not the reads from these libraries.

This distinction has to be explicit in the figure caption. A reader who sees a
tree in a 16S paper will assume the sequences came from the samples, and a
reviewer who discovers otherwise mid-review will treat it as a misrepresentation
rather than a design choice. The honest wording is along the lines of:
*"Reference 16S rRNA sequences (NCBI) for the N most abundant species detected,
aligned with MAFFT and trimmed with trimAl; maximum-likelihood tree inferred
with IQ-TREE under GTR+G with 1,000 ultrafast bootstrap replicates. Rings show
phylum, dominant substrate, and ESKAPEE membership."*

The tree therefore shows how the detected taxa are related to one another. It
does not show strain-level relationships within this dataset, and nothing in it
can support a claim about which lineages were transported along the transect.

## Reproducibility: the accession table

`01` writes `output/phylo_tip_provenance.csv`, one row per tip:

| column | meaning |
|---|---|
| `query_species` | the name the tip is labelled with |
| `accession` | the NCBI accession actually retrieved |
| `fetched_as` | `species` or `genus_fallback` |
| `retrieved_title` | the record's own description |

**Deposit this table.** The searches use `sort="relevance"`, so NCBI can return
a different top hit next year; without the accessions the tree cannot be
rebuilt, by a reviewer or by you.

## The genus fallback, and why the table matters

When no species-level record is found, the original notebook fell back to a
type-strain record **of the genus** and then relabelled it with the species
name. A tip reading *Legionella tunisiensis* could carry a different
*Legionella* sequence entirely, and nothing recorded that it had happened.

The fallback is kept — dropping those taxa would thin the tree — but every
fallback is now recorded in `fetched_as`, which the annotation step reads so
that those tips can be marked. Species-level resolution cannot be claimed for
them.

## Changes from the original notebook

1. **Accessions and fallback status are recorded.** They were not, which made
   the tree unreproducible and the fallback invisible.
2. **Plastic and water columns are selected by position** (`c[2]`), matching
   the rest of the pipeline. `"P" in c` matches any column containing the
   letter P anywhere and is fragile to renaming.
3. **trimAl failure is fatal** rather than silently falling through to the
   untrimmed alignment. The silent fallback would have put a Methods statement
   in the manuscript that did not describe what was run.

## To report in the Methods

Fill from the run:

```
sequences retrieved   N of M species (K by genus fallback)
mafft                 --auto  (version: CHECK)
trimal                -automated1  (version: CHECK)
iqtree                -m GTR+G -B 1000 -T AUTO  (version: CHECK)
alignment             CHECK columns before trimming, CHECK after
```

`iqtree` reports its version and the retained alignment length in the `.log`
and `.iqtree` files — take the numbers from there, not from memory.

Consider `-m MFP` (ModelFinder) instead of a fixed `GTR+G`: it costs little and
removes an obvious reviewer question. If GTR+G is kept, say why in the Methods.
