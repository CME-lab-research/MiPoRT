# MiPORT analysis workflow (STORM chart source)

Purpose:
Internal documentation of the MiPORT data processing and downstream analysis workflow used for workflow visualization and manuscript figure preparation.

Project:
MiPORT: Multi-cohort respiratory metagenome atlas

Last updated:
2026-05-06

Maintainers:
Tejus Shinde et al.

Notes:
- Sample counts reflect post-filtering datasets used in manuscript analyses.
- This document is intended for internal tracking and figure generation.
- Downstream analyses may evolve during manuscript revision.

```text
Data analysis workflow
│
├── Respiratory metagenome collection
│   ├── Public metagenomes (N = 4,123)
│   ├── Contributed metagenomes (N = 330)
│   └── Total metagenomes collected (N = 4,453)
│
├── Metadata harmonization and sample filtering
│   ├── Collect and standardize sample metadata
│   ├── Assign harmonized variables
│   │   ├── SampleType
│   │   ├── RT category
│   │   ├── Disease
│   │   ├── Healthy status
│   │   ├── BioProject / cohort
│   │   └── Other available covariates
│   └── Metadata-based filtering
│       └── Excluded samples (N = 976)
│
├── Uniform metagenome processing
│   ├── Read quality control
│   │   └── QC-failed metagenomes excluded (N = 203)
│   │
│   ├── Taxonomic profiling
│   │   └── Species-level profiles
│   │       └── Profiling-failed metagenomes excluded (N = 139)
│   │
│   ├── Functional profiling
│   │   ├── Gene family profiling
│   │   └── Regrouping to KEGG orthologs (KOs)
│   │
│   ├── Antibiotic resistance gene profiling
│   │   └── ARG abundance profiles
│   │
│   ├── Genome-resolved metagenomics
│   │   └── Metagenome assembly and MAG recovery
│   │
│   └── Batch-effect correction
│
└── Downstream analyses (N = 3,135 metagenomes)
    │
    ├── Taxonomic profiles (species-level)
    │   │
    │   ├── Alpha-diversity analysis
    │   │
    │   ├── Beta-diversity analysis
    │   │   ├── Ordination by SampleType and RT category
    │   │   ├── Within-SampleType β-diversity comparisons
    │   │   └── Variance partitioning by metadata variables (PERMANOVA/Adonis)
    │   │
    │   ├── Shared species analysis
    │   │   └── UpSet plot
    │   │
    │   ├── Differential prevalence and abundance testing
    │   │   └── MaAsLin3
    │   │       │
    │   │       ├── Lower respiratory tract
    │   │       │   └── BAL samples (N = 292; 5 datasets)
    │   │       │       ├── Healthy (N = 32) vs Pneumonia (N = 154)
    │   │       │       └── Healthy (N = 32) vs COVID-19 (N = 105)
    │   │       │
    │   │       └── Intermediate respiratory tract
    │   │           └── Sputum samples
    │   │               └── Healthy (N = 14) vs Cystic fibrosis (N = 138)
    │   │
    │   └── Microbial source contribution analysis
    │       └── FEAST
    │           ├── BAL sinks (N = 345) with RT sources (N = 956)
    │           └── Sputum sinks (N = 68) with RT sources (N = 960)
    │
    ├── Functional profiles (KEGG orthologs)
    │   │
    │   ├── Gene family regrouping
    │   │   └── Regrouped to KEGG orthologs (KOs)
    │   │
    │   ├── Beta-diversity analysis
    │   │   ├── Ordination by SampleType and RT category
    │   │   └── Variance partitioning by metadata variables (PERMANOVA/Adonis)
    │   │
    │   └── Differential prevalence and abundance testing
    │       └── MaAsLin3
    │           │
    │           ├── Lower respiratory tract
    │           │   └── BAL samples
    │           │       ├── Healthy (N = 32) vs Pneumonia (N = 154)
    │           │       └── Healthy (N = 32) vs COVID-19 (N = 105)
    │           │
    │           └── Intermediate respiratory tract
    │               └── Sputum samples
    │                   └── Healthy (N = 14) vs Cystic fibrosis (N = 138)
    │
    ├── Antibiotic resistance gene profiles
    │   │
    │   ├── ARG composition and clustering analysis
    │   │   └── Complex heatmap
    │   │
    │   └── Total ARG abundance across human RT sample types
    │       └── Bubble plot
    │
    └── Genome-resolved metagenomics
        │
        ├── Input metagenomes for assembly
        │   └── N = 3,135 samples
        │
        ├── Metagenome assembly and binning
        │   └── Recovered genomic bins (N = 12,609)
        │
        ├── Genome dereplication
        │   └── Dereplicated at 95% ANI
        │
        ├── RESPICAT genome catalog
        │   └── 1,255 non-redundant species-level genomes
        │
        ├── Taxonomic composition of recovered genomes
        │   └── Treemap by microbial order and phylum
        │
        └── Genome quality and assembly metrics
```

## Sample accounting

Total collected metagenomes: 4,453
- Metadata-filtered exclusions: 976
- QC-failed exclusions: 203
- Additional downstream exclusions: 139

Final downstream analysis set: 3,135
