# MiPoRT
MiPoRT (Microbial Profiling of the Respiratory Tract) study aims to study the microbial profiles and their ecology in the human respiratory tract. 

# Folder structure
│
├── README.md
├── requirements.text			# contains list of tools and package used with their versions
├── LICENSE
├── .gitignore			# specify which files not to upload
│
├── data/
│   ├── raw/                 # Feature tables from microbial profiles
│   ├── processed/           # Preprocessed and filtered sample profiles
│   └── metadata/            # Sample metadata
│
├── Analysis/
│   ├── 01_preprocessing.ipynb
│   ├── 02_alpha_beta_diversity.ipynb
│   ├── 03_differential_abundance.ipynb
│   └── 04_strain_analysis.ipynb
│
├── scripts/
│   ├── utils.py             # Common helper functions
│   ├── run_analysis.sh      # Bash scripts to execute pipelines
│   └── stats_model.R        # R script for regression or stat tests
│
├── results/
│   ├── tables/
│   └── figures/
│
└── docs/                    # Manuscript, supplementary text, notes