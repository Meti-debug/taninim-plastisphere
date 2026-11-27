# Input data

The scripts read their input tables from this folder. Only `metadata.tsv` is
kept in the repository; the count tables are deposited with the manuscript as
Supplementary Data. Raw reads are in the NCBI SRA under BioProject
PRJNA1458552.

| file | contents | used by |
|---|---|---|
| `counts.csv` | Kraken2/Bracken species counts from the EPI2ME `wf-metagenomics` workflow (`Lig_wf-metagenomics-counts-species.csv`): `species`, one column per library, then `kingdom` to `genus` and `tax` | 01, 02, 04, 05, 06 |
| `otu_table.tsv` | 97% OTU x library read counts, 24 libraries, reference-free (built from cluster consensus sequences) | 03 |
| `metadata.tsv` | one row per library: `SampleID` (barcode), `SampleName` (library label), `SampleType`, `Site`, `Environment`, `Replicate` | 02 (R scripts), 06 |

Place the two count tables here before running:

```
data/
├── counts.csv
├── otu_table.tsv
└── metadata.tsv
```

## Library labels

The count tables and `metadata.tsv` use the labels from sequencing. The
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

## Two taxonomic sources, never combined

- **Species-level profiles** (`counts.csv`): diversity, ordination,
  PERMANOVA, core taxa, pathogen abundances, source tracking and the
  reference tree.
- **97% OTUs** (`otu_table.tsv`): OTU richness and the shared and unique
  comparisons between sites and sample types only. No OTU carries a taxonomic
  assignment, and no Kraken2 species is treated as an OTU.
