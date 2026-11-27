# 04 Core taxa

Core species per habitat and sample type, and how many river-core species are
retained in the sea; plus the 16S reference tree of the detected species.

Generates:
- Figure 4D (core taxa) and the core-retention test
- Figure 4C (16S reference tree, drawn and annotated in iTOL)

Input:
- data/counts.csv
- NCBI (reference 16S sequences, fetched by the tree scripts)

Output:
- `script.py` prints the core lists and the permutation test to the console
- analysis/04_core_taxa/output/phylo_* (tree, alignment, tip provenance, tool versions)

Run from the repository root:

```bash
python analysis/04_core_taxa/script.py

export NCBI_EMAIL=you@example.org
bash analysis/04_core_taxa/reference_16S_tree/run_all.sh
```

Core definition: relative abundance > 0.1% in at least 5 of the 6 samples of a
habitat (river = R1 + R2, marine = S1 + Sc), separately for plastic-bound and
planktonic samples. Retention is the share of river-core species that are also
marine-core; plastic and water are compared by shuffling sample type within
site (9,999 shuffles). See the header of `script.py`.

The tree steps and the accession table are described in
`reference_16S_tree/README.md`.
