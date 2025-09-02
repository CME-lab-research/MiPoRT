# MiPoRT (Microbial Profiling of the Respiratory Tract)
MiPoRT study aims to study the microbial profiles and their ecology in the human respiratory tract. 

## Abstract
The MiPoRT (Microbial Profiling of the Respiratory Tract) project integrates over 4,000 metagenomic samples from the human respiratory tract to profile microbial communities across health and disease. This repository hosts the code, workflows, and processed outputs used in the MiPoRT study. Raw sequencing data are archived at ENA under accession PRJEB96845. Supplementary processed data and tables are available via [Add zenodo DOI].

## 📁 Folder Structure

```text
MiPoRT/
│
├── README.md
├── requirements.txt		# contains list of tools and package used with their versions. Not an exhaustive list.
├── LICENSE
├── .gitignore				# specify which files not to upload
│
├── Data/
│   ├── Raw/				# Raw feature tables of microbial profiles
│   ├── Processed/			# Pre-processed and filtered microbial profiles
│   └── Metadata/			# Sample metadata
│
└──  Scripts/
	├── 01_Preprocessing					# Scripts for preprocessing reads
	├── 02_Diversity_analysis				# Microbiome diversity analysis scripts
	├── 03_Abundance_Prevalence_analysis	# Scripts for plots and stats
	├── 04_Differential_modeling    # Scripts for source tracking analysis
	└── 05_Source_tracking			# Scripts for source tracking analysis

```

> Note: The scripts in *Scripts/* are developed iteratively and are not optimized for general re-use. We provide them for transparency and reproducibility. Only the essential scripts that reproduce key analyses are provided. 

## Data availability

Raw reads: ENA cite [PRJEB96845](https://www.ebi.ac.uk/ena/browser/view/PRJEB96845) .
Reused data: List of the BioProjects pulled from ENA/SRA is in Supplementary Table-1.

