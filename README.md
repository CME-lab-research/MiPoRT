# MiPoRT (Microbial Profiling of the Respiratory Tract)
MiPoRT study aims to study the microbial profiles and their ecology in the human respiratory tract. 

This repository contains code and files to be shared for the research paper titled: **'Pan-microbiome analysis along the human respiratory axis reveals an ecological continuum in health and disruption in disease'** 

## Abstract
The human respiratory tract (RT) harbors spatially structured microbial communities, yet their organization across anatomical sites and diseases remains poorly resolved. Here, we present a respiratory microbiome atlas integrating >4,000 metagenomes from upper, intermediate, and lower RT compartments across health, pneumonia, COVID-19, and cystic fibrosis. Standardized taxonomic and functional profiling revealed strong biogeographic organization, with healthy lower RT communities representing filtered subsets of upper RT microbiota. Disease was associated with reproducible depletion of core commensals, altered inter-compartmental microbial connectivity, and increased contributions from taxa of unknown origin. Functional analyses identified conserved repertoires across upper compartments and niche-specific enrichment in the lower RT. Prevalence-based models outperformed abundance-based approaches in detecting disease-associated microbiome disruption. Genome-resolved analyses recovered 1,255 dereplicated metagenome-assembled genomes, establishing a reusable respiratory genome catalog. Together, this atlas provides a framework for studying respiratory microbiome organization and dysbiosis across health and disease.

## 📁 Folder Structure

```text
MiPoRT/
│
├── README.md
├── requirements.txt                         # Contains a list of tools and packages used with their versions. Not an exhaustive list.
├── LICENSE
├── .gitignore                               # Specifies files not to upload
│
├── Data/
│   ├── Raw/                                 # Raw feature tables of microbial profiles
│   ├── Processed/                           # Pre-processed and filtered microbial profiles
│   └── Metadata/                            # Sample metadata
│
├── Docs/
│   ├── workflows/
│   │   └── MiPORT_analysis_workflow_STORMchart.md
│   │
│   └── BC_ordination_files/                 # Interactive Bray-Curtis PCoA HTML files
│
├── Manuscript_plots/                        # Raw image files from the manuscript
│   ├── Figure_1/
│   ├── Figure_2/
│   ├── Figure_3/
│   ├── Figure_4/
│   ├── Figure_5/
│   ├── Figure_6/
│   └── Supplementary/
│
├── Supplementary materials/
│   ├── MiPoRT Supplementary Materials.pdf   # Supplementary figures accompanying the manuscript
│   └── MiPoRT_Supplementary_Data-1.xlsx     # Supplementary tables accompanying the manuscript
│
└── Scripts/
    ├── 01_Preprocessing                     # Scripts for preprocessing reads
    ├── 02_Diversity_analysis                # Microbiome diversity analysis scripts
    ├── 03_Abundance_Prevalence_analysis     # Scripts for plots and stats
    ├── 04_Differential_modeling             # Scripts for differential modeling analysis
    └── 05_Source_tracking                   # Scripts for source tracking analysis
    └── 06_Functional_analysis/
        ├── 01_A_process_humann_outputs.sh                    # Merge, clean, and regroup HUMAnN-derived gene-family and pathway profiles
        ├── 01_B_filter_zero_KO_samples.R                     # Remove samples or KO features with zero functional signal before downstream normalization/modeling
        ├── 02_A_Batch_Correction_split_KO_table_by_RT-Cat.R   # Split KO abundance tables by respiratory tract category before within-compartment batch correction
        ├── 02_B_Batch_Correction_withinRT_MMUPHin.R           # Perform within-respiratory-tract batch correction of KO profiles using MMUPHin
        ├── 02_C_Merge_batch_correction_tables.R               # Merge batch-corrected KO tables from different respiratory tract compartments
        ├── 03_pcoA_KO_ext_RELAB.R                             # Generate PCoA ordinations from relative-abundance KO profiles
        ├── 04_Permanova_analysis_with_Adonis_KO_relab.R       # Test functional profile variation using PERMANOVA/adonis on KO relative-abundance data
        └── 05_Functional_BRITE_heatmap-bubble_plot_script.R   # Generate KEGG BRITE-level heatmap and bubble plots for functional summaries
        
```

> Note: The scripts in *Scripts/* are developed iteratively and are not optimized for general re-use. We provide them for transparency and reproducibility. Only the essential scripts that reproduce key analyses are provided. 

## Repository contents

This repository contains custom code, analysis workflows, manuscript plots, supplementary files, and processed microbial profiles associated with the MiPoRT respiratory metagenome atlas.

- `Data/` contains raw and processed microbial feature tables, together with sample metadata used for downstream analyses.
- `Docs/` contains workflow documentation and interactive ordination files.
- `Manuscript_plots/` contains raw image files for the main manuscript figures and supplementary figures.
- `Supplementary materials/` contains the supplementary figures PDF and supplementary data table file accompanying the manuscript.
- `Scripts/` contains scripts used for preprocessing, diversity analysis, abundance/prevalence analysis, differential modeling, and source-tracking analyses.

## Data availability

The raw sequencing reads generated and analyzed in this study have been deposited in the European Nucleotide Archive (ENA) under project accession [PRJEB96845](https://www.ebi.ac.uk/ena/browser/view/PRJEB96845).

Previously published datasets were re-analysed from the ENA/NCBI Sequence Read Archive. These include BioProject accessions `PRJNA48479` and `PRJNA917836`, among others. A complete list of all reused BioProjects, including sample counts and sequencing statistics, is provided in Supplementary Table 2 in `MiPoRT_Supplementary_Data-1.xlsx`.

RESPICAT, the non-redundant metagenome-assembled genome catalogue generated in this study, is archived on Zenodo at [https://doi.org/10.5281/zenodo.20162797](https://doi.org/10.5281/zenodo.20162797). Accompanying documentation, checksums, and download/verification scripts are available through the RESPICAT GitHub repository: [https://github.com/CME-lab-research/RESPICAT](https://github.com/CME-lab-research/RESPICAT).

Custom code and analysis workflows for the MiPoRT atlas are available from this repository: [https://github.com/CME-lab-research/MiPoRT](https://github.com/CME-lab-research/MiPoRT).

Taxonomic, functional, and antimicrobial-resistance-gene (ARG) profiles generated for the MiPoRT atlas are also available through this repository (See [`Data/Processed/`](Data/Processed/) directory).

## Citation

Citation information will be added after publication or public release of the corresponding manuscript/preprint.
