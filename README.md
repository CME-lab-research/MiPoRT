# MiPoRT
MiPoRT (Microbial Profiling of the Respiratory Tract) study aims to study the microbial profiles and their ecology in the human respiratory tract. 

# 📁 Folder Structure

```text
MiPoRT/
│
├── README.md
├── requirements.text			# contains list of tools and package used with their versions
├── LICENSE
├── .gitignore			# specify which files not to upload
│
├── Data/
│   ├── raw/			# Feature tables from microbial profiles
│   ├── processed/			# Preprocessed and filtered sample profiles
│   └── metadata/			# Sample metadata
│
├── Analysis/
│   ├── 01_preprocessing
│   ├── 02_batchCorrection
│   ├── 03_alpha_beta_diversity
│   └── 04_Differential_analysis
│
├── Scripts/
│   ├── Data_download			# Common helper functions
│   ├── QC_preproc			# Common helper functions
│   ├── run_analysis.sh			# Bash scripts to execute pipelines
│   └── stats_model.R			# R script for regression or stat tests
│
├── Results/
│   ├── tables/
│   └── figures/
│
└── docs/			# Manuscript, supplementary text, notes
```

