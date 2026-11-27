# The Taninim River-Sea Continuum: Plastisphere’s Persistence, Connectivity, and Potential Pathogen Dispersal into the Eastern Mediterranean Sea

Analysis code.

- **Species-level profiles** from Kraken2 and Bracken, run through the EPI2ME
  `wf-metagenomics` workflow. Taxonomy assignment comes
  from here: diversity, ordination, PERMANOVA, core taxa, pathogen abundances,
  source tracking and the reference tree.
- **97% OTUs**, reference-free, built from cluster consensus sequences. Used
  only for OTU richness and for the shared and unique comparisons between
  sites and sample types. No OTU carries a taxonomic assignment and no Kraken2
  species is treated as an OTU; the two are never combined within a single
  statistic.

## Layout

```
taninim-plastisphere/
├── README.md
├── requirements.txt          Python packages
├── environment.yml           conda environment (Python, R and command-line tools)
├── data/                     input tables; see data/README.md
├── figures/                  final Figures 1-5 and supplementary figures
├── analysis/
│   ├── paths.R               shared paths for the R scripts
│   ├── 01_alpha_diversity/   Shannon diversity                      Fig. 3A-E
│   ├── 02_beta_diversity/    Bray-Curtis, PCoA, PERMANOVA           Fig. 3F-H, Table S1
│   ├── 03_OTU/               97% OTU persistence and sensitivity    Fig. 4A-B
│   ├── 04_core_taxa/         core taxa and the 16S reference tree   Fig. 4C-D
│   ├── 05_pathogens/         pathogen and ESKAPEE abundance         Fig. 5A-B
│   ├── 06_FEAST/             source tracking and sensitivity        Fig. 4E, Table S6, Fig. S7
│   └── 07_CTXM/              blaCTX-M group 1 tree                  Fig. 5C, Fig. S6
└── supplementary/            supplementary tables
```

Each analysis folder has a `README.md` listing what it generates, what it
reads and what it writes, and an `output/` folder for its results. Outputs are
not tracked in git, except the deposited CTX-M tree in `analysis/07_CTXM/tree/`.

## Running

Put `counts.csv` and `otu_table.tsv` in `data/` (see `data/README.md`), then,
from the repository root:

```bash
conda env create -f environment.yml
conda activate taninim-plastisphere
Rscript -e 'remotes::install_github("cozygene/FEAST")'

python  analysis/01_alpha_diversity/script.py
python  analysis/02_beta_diversity/script.py
Rscript analysis/02_beta_diversity/dispersion_turnover.R
Rscript analysis/02_beta_diversity/metric_robustness.R
python  analysis/03_OTU/script.py data/otu_table.tsv
python  analysis/04_core_taxa/script.py
python  analysis/05_pathogens/script.py
Rscript analysis/06_FEAST/primary_analysis.R
Rscript analysis/06_FEAST/sensitivity_analysis.R
Rscript analysis/06_FEAST/genus_level.R
Rscript analysis/06_FEAST/river_signature.R
Rscript analysis/06_FEAST/marine_site_comparison.R
Rscript analysis/07_CTXM/02_plot_complete_tree.R
Rscript analysis/07_CTXM/03_plot_figure5c.R
Rscript analysis/07_CTXM/04_plot_figure_s6.R
```

Panel C as placed in the composite Figure 5:

```bash
CTXM_STEM=Figure5C_CTXM_composite CTXM_W=2.93 CTXM_H=1.74 CTXM_TEXT=0.66 \
CTXM_TREE_FR=0.40 CTXM_ROOT_LEN=0.013 Rscript analysis/07_CTXM/03_plot_figure5c.R
```

Two steps need data that are not in the repository: the 16S reference tree
(`analysis/04_core_taxa/reference_16S_tree/run_all.sh`) fetches sequences from
NCBI, and rebuilding the CTX-M tree
(`analysis/07_CTXM/01_build_panel_and_tree.sh`) needs the amplicon consensus
and reference sequences in `$CTXM_DATA`.

## Requirements

Python 3.11 or later with numpy, pandas, scipy and biopython
(`requirements.txt`). R 4.4 or later with vegan, ggplot2, dplyr, tidyr,
readr, tibble, scales, ape, ggtree, treeio and FEAST
(github.com/cozygene/FEAST, v0.1.0, not the Bioconductor package of the same
name). Command-line tools for the trees: MAFFT, trimAl, IQ-TREE 2 and seqkit;
the versions used for the CTX-M tree are in
`analysis/07_CTXM/tree/ctxm_versions.txt`. `environment.yml` installs all of
these except FEAST.
