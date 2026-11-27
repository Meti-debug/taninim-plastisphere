# Phylogeny of the detected taxa

data/counts.csv -> species list -> NCBI 16S -> MAFFT -> trimAl -> IQ-TREE

Everything at once, from the repository root:

```bash
export NCBI_EMAIL=metigizawmaranata@gmail.com
bash analysis/04_core_taxa/reference_16S_tree/run_all.sh
```

Or one step at a time:

| step | script |
|---|---|
| build the species list and fetch one 16S sequence per species | `01_fetch_reference_16S.py` |
| align, trim, infer the maximum-likelihood tree | `02_align_trim_tree.sh` |

The tree was drawn and annotated in iTOL v6 (branch colours by site, inner ring
by sample type, outer ring by family).

`trim_alignment.py` is a helper: it runs the `trimal` binary when it is on PATH
and pytrimal otherwise, because trimal is not packaged on Ubuntu 24.04.

Dependencies:

```bash
sudo apt-get install -y mafft iqtree
pip install biopython pandas numpy pytrimal
```

## What this tree is

The tips are **reference 16S sequences downloaded from NCBI for each species
name detected in this study**. They are not the reads from these libraries.

Reference 16S rRNA sequences (NCBI) for the N most abundant species detected,
aligned with MAFFT and trimmed with trimAl; maximum-likelihood tree inferred
with IQ-TREE under GTR+G with 1,000 ultrafast bootstrap replicates. Rings show
the most dominant family.

## Reproducibility: the accession table

`01` writes `analysis/04_core_taxa/output/phylo_tip_provenance.csv`, one row per tip:

| column | meaning |
|---|---|
| `query_species` | the name the tip is labelled with |
| `accession` | the NCBI accession actually retrieved |
| `fetched_as` | `species` or `genus_fallback` |
| `retrieved_title` | the record's own description |

without the accessions the tree cannot be
reproduced.

## Methods

Fill from the run:

```
sequences retrieved   N of M species (K by genus fallback)
mafft                 --auto  (version: CHECK)
trimal                -automated1  (version: CHECK)
iqtree                -m GTR+G -B 1000 -T AUTO  (version: CHECK)
alignment             CHECK columns before trimming, CHECK after
