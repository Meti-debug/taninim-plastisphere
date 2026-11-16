# The Taninim River-Sea Continuum: Plastisphere’s Persistence, Connectivity, and Potential Pathogen Dispersal into the Eastern Mediterranean Sea

Analysis code.

- **Species-level profiles** from Kraken2 and Bracken, run through the EPI2ME
  `wf-metagenomics` workflow. Taxonomy assignment comes
  from here: diversity, ordination, PERMANOVA, core taxa, pathogen abundances,
  source tracking and the reference tree.
- **97% OTUs**, reference-free, built from cluster consensus sequences. 

## Requirements

```
python >= 3.11        R >= 4.4
numpy >= 2.0          ggplot2 >= 3.5
pandas >= 2.0         pheatmap >= 1.0.12
scipy >= 1.11         dplyr, tidyr, scales
biopython >= 1.83     (phylogeny only)
```

External tools, for `analysis/phylogeny/` only: MAFFT v7, trimAl, IQ-TREE 2 v2.3.5.


