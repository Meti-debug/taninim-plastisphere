# Input data

The scripts read their input tables from this folder. The tables are not kept
in the repository; they are deposited with the manuscript as Supplementary
Data. Raw reads are in the NCBI SRA under BioProject
PRJNA1458552.

| file | contents | used by |
|---|---|---|
| `counts.csv` | Kraken2/Bracken species counts from the EPI2ME `wf-metagenomics` workflow (`Lig_wf-metagenomics-counts-species.csv`): `species`, one column per library named by its library label, then the taxonomy columns | 01, 02, 04, 05, 06 |
| `otu_table.tsv` | 97% OTU x library read counts, 24 libraries, reference-free (built from cluster consensus sequences) | 03 |

Place the two count tables here before running:

```
data/
├── counts.csv
└── otu_table.tsv
```

Site and sample type are read from the library labels in the column names, so
no separate sample sheet is needed.

## Library labels

The count tables use the labels from sequencing. The
scripts convert them to the manuscript labels:

| sequencing | manuscript |
|---|---|
| R2 | R1 (Ada confluence) |
| R3 | R2 (estuary) |
| S1P4 | S1P1 |
| S1, Sc | unchanged |

Library names are site + sample type + replicate, e.g. `R3P2` is the second
plastic-bound replicate from the estuary. P = plastic-bound, W = water
(planktonic).

The EPI2ME export names its columns by barcode. The barcodes correspond to
these library labels:

| barcode | 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| library | R2W1 | R2W2 | R2W3 | R2P1 | R2P2 | R2P3 | R3W1 | R3W2 | R3W3 | R3P1 | R3P2 | R3P3 |

| barcode | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| library | S1W1 | S1W2 | S1W3 | S1P4 | S1P2 | S1P3 | ScW1 | ScW2 | ScW3 | ScP1 | ScP2 | ScP3 |

## Two taxonomic sources, never combined

- **Species-level profiles** (`counts.csv`): diversity, ordination,
  PERMANOVA, core taxa, pathogen abundances, source tracking and the
  reference tree.
- **97% OTUs** (`otu_table.tsv`): OTU richness and the shared and unique
  comparisons between sites and sample types only. No OTU carries a taxonomic
  assignment, and no Kraken2 species is treated as an OTU.
