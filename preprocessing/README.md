# Preprocessing: reads to abundance tables

Two tables come out of this folder, and the rest of the repository starts from
them.

    pass reads -> EPI2ME (Kraken2 + Bracken) -> data/counts.csv
               -> QC filter -> consensus clustering -> otu_table.tsv

| step | script |
|---|---|
| taxonomic classification | `01_epi2me_kraken2.md` (a record, not a script) |
| EPI2ME counts to a tidy table | `02_epi2me_counts_to_table.py` |
| 97% OTU table from the same reads | `03_otu97_pipeline.sh` |

Basecalling and demultiplexing ran live in MinKNOW on the sequencer; the
instrument, kit, basecaller version and read counts are in the manuscript
Methods, and the pass reads are deposited in the SRA, so the analysis starts
from those reads rather than from the raw signal.

Classification was run in the EPI2ME desktop application rather than from a
shell, so `01_epi2me_kraken2.md` records the workflow version, the database and
every parameter instead of pretending to be a script.

`03_otu97_pipeline.sh` is self-contained and was run from a shell. It writes
`run_settings.txt` recording the VSEARCH version and every parameter used, and
`read_tracking.tsv` recording reads retained at each step, so a run can be
audited after the fact.
