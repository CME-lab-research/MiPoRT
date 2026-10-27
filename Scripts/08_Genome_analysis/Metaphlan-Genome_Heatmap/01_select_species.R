# Code to select species to plot figure: Respiratory niche distribution and disease-linked functional potential of RESPICAT-represented taxa.

# load libs
library(tidyverse)

# 1. Read RESPICAT genome taxonomy file
RESPICAT_genome_taxonomy <- read_tsv(file = "./01_Select_species/RESPICAT_MAG_taxonomy_v1.0.tsv")

    # check how many have species assignments
    RESPICAT_genome_taxonomy %>%
        count(no_species_assignment)
    
    "
    # A tibble: 2 × 2
      no_species_assignment     n
      <lgl>                 <int>
    1 FALSE                   651
    2 TRUE                    604
    "

    # Select genomes with species assignments
    RESPICAT_genome_taxonomy <- 
        RESPICAT_genome_taxonomy %>%
        filter(no_species_assignment == FALSE)

    # Retained 651 genomes which have species annotation
    
# 2. Read Maaslin-3 full model files for LRT and IRT
tmp_LRT_species <- read_tsv("01_Select_species/LRT_plot_df.txt")
tmp_IRT_species <- read_tsv("01_Select_species/IRT_plot_df.txt")

# Add sample type information as a column
tmp_LRT_species$SampleType = "BAL"
tmp_IRT_species$SampleType = "Sputum"

    # Merge 2 tables
    # select cols found in both tables
    common_cols <- intersect(
    colnames(tmp_LRT_species),
    colnames(tmp_IRT_species))
    
    # row bind 
    Disease_associated_species <- rbind(
        tmp_LRT_species[,common_cols],
        tmp_IRT_species[,common_cols])
    
    # mv SampleType column to front
    Disease_associated_species <- Disease_associated_species %>%
        filter(ID == "Model1", metadata == "Disease") %>%
        relocate(all_of(c("model", "SampleType"))) %>%
        group_by(SampleType)
    
    # filter for high confidence association with p-value and q-value cutoffs. q < 0.05 and |β| ≥ 1
    Disease_associated_species_filt <- 
        Disease_associated_species %>%
        filter(qval_individual < 0.05)

    Disease_associated_species_filt %>%
        count(SampleType, model, value)

    # Clean species name column
    # Create a new column for species names
    
    # Load the helper function
    source("C:/Users/o_shinde/Desktop/MiPORT_DA/Scripts/R_Utility_scripts/format_species_label.R")  # adjust path if needed
    
    # Format variable 'feature' and  create a new col 'Species names'
    Disease_associated_species_filt$Species <- 
        format_species_label(
            Disease_associated_species_filt$feature,
            italicize = FALSE)
    
    # relocate species names to front
    Disease_associated_species_filt <- 
        Disease_associated_species_filt %>% 
        relocate(Species, .before = feature) %>% 
        relocate(feature, .after = last_col())
    
    # count of unique species with any association
    unique(Disease_associated_species_filt$Species)
        # 102

# 3. Check intersection of RESPICAT and Maaslin-3 species discovery
    which(Disease_associated_species_filt$Species %in%
              RESPICAT_genome_taxonomy$species)
    
    # names of common species
    Species_to_Retain <- intersect(
        Disease_associated_species_filt$Species,
        RESPICAT_genome_taxonomy$species
    )
    
    # sort the list
    Species_to_Retain <- sort(Species_to_Retain)
    
    # subset Maaslin-3 results for these species
    Select_disease_associated_species <- 
        Disease_associated_species_filt %>%
        filter(Species %in% Species_to_Retain)
    
    #save species names to highlight
    write_tsv(Select_disease_associated_species,
              "Maaslin-3-associated_RESPICAT_intersected_Species.txt")
    
    
    
    
    
    
    
    
    