# Input tables

Place these in `data/` before running anything.

| file | contents | written by |
|---|---|---|
| `counts.csv` | Kraken2/Bracken species counts: `tax`, `genus`, `species`, then one column per library | `preprocessing/02_epi2me_counts_to_table.py` |
| `otu_table.tsv` | 97% OTU x library read counts, 24 libraries | `preprocessing/03_otu97_pipeline.sh` |

Both are deposited with the manuscript as Supplementary Data. Raw reads are in
the NCBI SRA under BioProject PRJNA1458552.

## Library labels

`counts.csv` carries the site labels used during sequencing, which are one
position upstream of the labels used in the manuscript:

| in `counts.csv` | in the manuscript |
|---|---|
| `R2P1`–`R2W3` | `R1P1`–`R1W3` |
| `R3P1`–`R3W3` | `R2P1`–`R2W3` |
| `S1P4` | `S1P1` |
| all others | unchanged |

The table is deposited with the sequencing labels, and the analysis scripts
convert them on read. The OTU table is written with manuscript labels already.
