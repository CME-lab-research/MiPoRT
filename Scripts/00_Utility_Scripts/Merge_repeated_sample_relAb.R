# Script to account for repeated measures.
# This script reads metadata and the metaphlan counts file, then finds samples that needs to be combined by taking a mean of relative abundance for all the bugs

library(tidyverse)

# Init vars
Data_dir="C:/Users/o_shinde/Desktop/MiPORT/AG Moissl-Eichinger_people/Tejus/MiPoRT Paper/MiPORT_DA/"

# Read metadata file
Metadata_Df <- read.table("./MiPORT_Metadata_downstream_filtered_M4_passed.tsv", header = T, sep = "\t")

# Make a list of metadata cols to work with
Groups_Interest <- c("PatientID", "SampleID", "BioProject", "Healthy", "SampleTypev2", "RT_category", "Disease")

# Read feature table
Metaphlan_profiles <- read.table(paste0(Data_dir, "MetaPhlan4_results/Metaphlan4_species_merged_filtered_batchCorrected.tsv"), header = T, sep = "\t")

# Read file to filter species present in min 10 samples
Species_filter_df <- read.table(paste0(Data_dir,"MetaPhlan4_results/Species_min_10samples.txt"), header = T, sep = "\t")

# Part a. Healthy

    # subset for healthy
    Metadata_Df_Healthy <- Metadata_Df %>%
        select(all_of(Groups_Interest)) %>%
        filter(Healthy == TRUE)
    
    Unique_individual_IDs <- unique(Metadata_Df_Healthy$PatientID)    
    # We have 754 unique individuals
    
    # Find samples present in both files to retain them
    Retain_Samples <- intersect(colnames(Metaphlan_profiles), Metadata_Df_Healthy$SampleID)
    
    # Filter metaphlan count table for samples of interest
    Metaphlan_profiles_Healthy <- Metaphlan_profiles %>% 
        rownames_to_column("taxonomy") %>%
        select(all_of(c("taxonomy", Retain_Samples))) %>% # filter samples
        filter(taxonomy %in% Species_filter_df$Species)
    
    # Filter metadata table for the same
    Metadata_Df_Healthy <- Metadata_Df_Healthy %>%
        filter(SampleID %in% Retain_Samples)
        # 1429 Healthy samples now, after correction there should be 
    
        table(Metadata_Df$SampleTypev2, useNA = 'ifany')
    
        # Step 1: Find patients with repeated samples *within each SampleType*
        Repeated_Sample_details <- Metadata_Df_Healthy %>%
            filter(!is.na(PatientID)) %>%
            group_by(SampleTypev2, PatientID) %>%
            summarise(Count_Patients = n(), .groups = "drop") %>%
            filter(Count_Patients > 1)
        
        #263 combinations
        
        # Step 2: Join to get all sample IDs from those SampleType-Patient combinations
        Samples_to_Merge <- Metadata_Df_Healthy %>%
            filter(!is.na(PatientID)) %>%
            inner_join(Repeated_Sample_details, by = c("SampleTypev2", "PatientID")) %>%
            arrange(SampleTypev2, PatientID)
        
        'sum(Repeated_Sample_details$Count_Patients)
    [1] 715'
    
        # 715 samples will be removed and merged into 263 samples.
    
    # Let's start merging
    # 1. Transpose your Metaphlan data to long format
    Metaphlan_long <- Metaphlan_profiles_Healthy %>%
        select(all_of(c("taxonomy", Samples_to_Merge$SampleID))) %>%
        pivot_longer(-taxonomy, names_to = "SampleID", values_to = "Abundance")
    
    # 2. Add PatientID using join
    Metaphlan_long_merged <- Metaphlan_long %>%
        left_join(Samples_to_Merge %>% select(SampleID, PatientID, SampleTypev2, Disease, RT_category), by = "SampleID")
    
    # 3. Calculate mean abundance per PatientID for each species
    Patient_mean_profiles <- Metaphlan_long_merged %>%
        group_by(PatientID, SampleTypev2, taxonomy) %>%
        summarise(Mean_Abundance = mean(Abundance), .groups = "drop") %>%
        mutate(NewSampleNames = paste(PatientID, SampleTypev2, sep = "_"))
    
    # 4. Reshape back to wide format (species x patient)
    Metaphlan_collapsed <- Patient_mean_profiles %>%
        select(NewSampleNames, taxonomy, Mean_Abundance) %>%
        pivot_wider(names_from = NewSampleNames, values_from = Mean_Abundance)
    
        # sanity check
        
    # 5. Remove the samples from the main df and replace them with these
    Corrected_Metaphlan_profiles_Healthy <- Metaphlan_profiles_Healthy %>%
        select(!all_of(Samples_to_Merge$SampleID)) %>% # Take samples not in Samples_to_Merge list. It has list of samples which were changed.
        cbind(Metaphlan_collapsed) %>% # This binds OG profiles to collapsed profiles (for samples in Samples_to_Merge list)
        column_to_rownames("taxonomy") %>%
        select(-taxonomy)
    
    '655 (samples w/o repeated measures) +
        715 samples to be merged into 322 samples, to give a total of 978 samples
    (655 + 160 samples to work with)'
    
    
    # save table
    write.table(Corrected_Metaphlan_profiles_Healthy,
                "../../MetaPhlan4_results/Repeated_Measure_merged/Healthy_Metaphlan4_species_merged_RepeatedSamples_batchCorrected.tsv",
                sep = "\t", row.names = T, quote = F)
    
    
    # save a metadata file for the same
    
    # 1. Get metadata for patients with repeated samples
    Repeated_metadata <- Samples_to_Merge %>%
        group_by(SampleTypev2, PatientID) %>%
        summarise(
            SampleID = first(SampleID),  # Keep the first sample ID (for reference, optional)
            BioProject = first(BioProject),  # Assuming same project for all samples
            Healthy = first(Healthy),        # Assuming same health status across replicates
            SampleTypev2 = first(SampleTypev2),  # Choose a strategy if they vary
            RT_category = first(RT_category)     # Same as above
        )
    
    # All patients
    Final_Patient_Metadata <- Metadata_Df_Healthy %>%
        filter(SampleID %in% Samples_to_Merge$SampleID) %>%
        group_by(SampleTypev2, PatientID) %>%
        summarise(
            SampleID = unique(paste(PatientID, SampleTypev2, sep = "_")),
            BioProject = first(BioProject),
            Healthy = first(Healthy),
            SampleTypev2 = first(SampleTypev2),
            RT_category = first(RT_category)
        ) %>%
        rbind(Metadata_Df_Healthy %>%
                  filter(!(SampleID %in% Samples_to_Merge$SampleID)) )
    
    write.table(Final_Patient_Metadata, "../../MetaPhlan4_results/Repeated_Measure_merged/Healthy_Metaphlan4_metadata_RepeatedSamples.tsv", sep = "\t", row.names = F, quote = F)
    
    write.table(Samples_to_Merge, "sample_Merge_index_RepeatedSamples.txt",
                row.names = F, quote = F, sep = '\t')

    
# Part b. Diseased
    # We take all samples > Healthy URT/IRT/LRT as source & Diseased LRT/IRT > Account for repeated measures > Save the table
    
    # Subset the main table as required for FEAST analysis
    
    # subset for Diseased samples only
    Metadata_Df_M <- Metadata_Df %>%
        select(all_of(Groups_Interest)) 
    
    Unique_individual_IDs <- unique(Metadata_Df_M$PatientID)    
    # We have 2049 unique individuals out of 2967 samples (2812 in both)
    
    # Find samples present in both files to retain them
    Retain_Samples <- intersect(colnames(Metaphlan_profiles), Metadata_Df$SampleID) # 2812
    
    # Filter metaphlan count table for samples of interest
    Metaphlan_profiles_D <- Metaphlan_profiles %>% 
        rownames_to_column("taxonomy") %>%
        select(all_of(c("taxonomy", Retain_Samples))) %>% # filter samples
        filter(taxonomy %in% Species_filter_df$Species)
    
    # Filter metadata table for the same
    Metadata_Df_M <- Metadata_Df_M %>%
        filter(SampleID %in% Retain_Samples)
    # 2812 samples in total: 1373 Diseased samples now, after correction there should be 1013?
    
    table(Metadata_Df_M$SampleTypev2, useNA = 'ifany')
    table(Metadata_Df_M$Healthy, useNA = 'ifany')
    
    # Step 1: Find patients with repeated samples *within each SampleType*
    Repeated_Sample_details <- Metadata_Df_M %>%
        filter(!is.na(PatientID)) %>%
        group_by(SampleTypev2, PatientID) %>%
        summarise(Count_Patients = n(), .groups = "drop") %>%
        filter(Count_Patients > 1)
    
    #376 combinations
    
    # Step 2: Join to get all sample IDs from those SampleType-Patient combinations
    Samples_to_Merge <- Metadata_Df_M %>%
        filter(!is.na(PatientID)) %>%
        inner_join(Repeated_Sample_details, by = c("SampleTypev2", "PatientID")) %>%
        arrange(SampleTypev2, PatientID)
    
    # 270 unique combinations from 1024 rows.
    
    length(unique(Repeated_Sample_details$PatientID))
    '[1] 273'
    
    # 1024 samples will be removed and merged into 273 samples.
    
    # Let's start merging
    # 1. Transpose your Metaphlan data to long format
    Metaphlan_long <- Metaphlan_profiles_D %>%
        select(all_of(c("taxonomy", Samples_to_Merge$SampleID))) %>%
        pivot_longer(-taxonomy, names_to = "SampleID", values_to = "Abundance")
    
    # 2. Add PatientID using join
    Metaphlan_long_merged <- Metaphlan_long %>%
        left_join(Samples_to_Merge %>% select(SampleID, PatientID, SampleTypev2, Disease, RT_category), by = "SampleID")
    
    # 3. Calculate mean abundance per PatientID for each species
    Patient_mean_profiles <- Metaphlan_long_merged %>%
        group_by(PatientID, SampleTypev2, taxonomy) %>%
        summarise(Mean_Abundance = mean(Abundance), .groups = "drop",
                  Disease = first(Disease),
                  RT_category = first(RT_category)) %>%
        mutate(NewSampleNames = paste(PatientID, SampleTypev2, sep = "_"))
    
    # 4. Reshape back to wide format (species x patient)
    Metaphlan_collapsed <- Patient_mean_profiles %>%
        select(NewSampleNames, taxonomy, Mean_Abundance) %>%
        pivot_wider(names_from = NewSampleNames, values_from = Mean_Abundance)
    
    # sanity check
    # 5. Remove the samples from the main df and replace them with these
    Corrected_Metaphlan_profiles_D <- Metaphlan_profiles_D %>%
        select(!all_of(Samples_to_Merge$SampleID)) %>%
        cbind(Metaphlan_collapsed) %>%
        column_to_rownames("taxonomy") %>%
        select(-taxonomy)
    
    '
    Final count of samples with repeated measures : 2164
    1891 samples are without repeated sampling.
    1024 samples to be merged into 273 samples, to give a total of 2164 samples : 1891 + 273 samples to work with.'
    
    # save table
    write.table(Corrected_Metaphlan_profiles_D,
                "./MetaPhlan4_results/All_Metaphlan4_species_merged_RepeatedSamples_batchCorrected.tsv",
                sep = "\t", row.names = T, quote = F)
    
    
    # save a metadata file for the same
    
    # 1. Get metadata for patients with repeated samples
    Repeated_metadata <- Samples_to_Merge %>%
        group_by(SampleTypev2, PatientID) %>%
        summarise(
            SampleID = first(SampleID),  # Keep the first sample ID (for reference, optional)
            BioProject = first(BioProject),  # Assuming same project for all samples
            Healthy = first(Healthy),        # Assuming same health status across replicates
            SampleTypev2 = first(SampleTypev2),  # Choose a strategy if they vary
            RT_category = first(RT_category),     # Same as above
            Disease = first(Disease)     # Same as above
        )
    
    # All patients
    Final_Patient_Metadata <- Metadata_Df_M %>%
        filter(SampleID %in% Samples_to_Merge$SampleID) %>%
        group_by(SampleTypev2, PatientID) %>%
        summarise(
            SampleID = unique(paste(PatientID, SampleTypev2, sep = "_")),
            BioProject = first(BioProject),
            Healthy = first(Healthy),
            SampleTypev2 = first(SampleTypev2),
            RT_category = first(RT_category),
            Disease = first(Disease)
        ) %>%
        rbind(Metadata_Df_M %>%
                  filter(!(SampleID %in% Samples_to_Merge$SampleID)) )
    
    write.table(Final_Patient_Metadata, "./MetaPhlan4_results/Repeated_Measure_merged/All_Metaphlan4_metadata_RepeatedSamples.tsv", sep = "\t", row.names = F, quote = F)
    
    write.table(Samples_to_Merge, "./MetaPhlan4_results/Repeated_Measure_merged/Sample_Merge_index_RepeatedSamples_diseased.txt",
                row.names = F, quote = F, sep = '\t')
    


    
    
