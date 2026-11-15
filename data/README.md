# Input tables

Place these in `data/` before running anything.

| file | contents |
|---|---|
| `counts.csv` | Kraken2/Bracken species counts from EPI2ME `wf-metagenomics`: `tax`, `genus`, `species`, then one column per library |
| `otu_table.tsv` | 97% OTU x library read counts, 24 libraries |

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
