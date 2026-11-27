# The Taninim River-Sea Continuum: Plastisphere’s Persistence, Connectivity, and Potential Pathogen Dispersal into the Eastern Mediterranean Sea

Analysis code.


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

