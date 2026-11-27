# 03 OTU

97% OTUs, used only for richness and for shared and unique OTUs between sites
and sample types. No OTU carries a taxonomic assignment.

Generates:
- the four-site persistence test (OTUs present at all four sites, plastic vs planktonic)
- the discovery-depth sensitivity comparison (supplementary)

Input:
- data/otu_table.tsv (OTU x library read counts, 24 libraries)

Output:
- both scripts print their results to the console

Run from the repository root:

```bash
python analysis/03_OTU/script.py      data/otu_table.tsv
python analysis/03_OTU/sensitivity.py TABLE_A.tsv TABLE_B.tsv [LABEL_A LABEL_B]
```

| script | what it does |
|---|---|
| `script.py` | share of OTUs present at all four sites, plastic minus water, within-site permutation test (9,999 shuffles) |
| `sensitivity.py` | compares two OTU tables built at different discovery depths: totals, site-specific, all-four-site, plastic-only / shared / water-only, richness |

Presence rule (as in Figure 4A-B): an OTU is present at a site, in one sample
type, when it has at least 8 reads in at least 2 of the 3 replicates.
