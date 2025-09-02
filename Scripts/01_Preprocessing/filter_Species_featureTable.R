# Script to filter merged Metaphlan feature table

# Load libraries
library(tidyverse)

# Read feature table & metadata table
Metaphlan_profiles <- read.table("../01_MMuphin_Batch_Correction/MiPORT_All_RT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", header = T, sep = "\t")
Metadata_Df <- read.table("../MiPoRT_Metadata_global_QC_passed_M4_profiled.tsv", header = T, sep = "\t")

# Identify samples with 0 total relative abundance
SampleFilter_colSums <- colSums(Metaphlan_profiles)
    # plot and see
    hist(SampleFilter_colSums, breaks = 50, main = "Sample total relative abundances",
     xlab = "Total abundance", col = "skyblue")
    abline(v = 0.001, col = "red", lty = 2)

    # sanity check
    length(which(SampleFilter_colSums == 0))
        # 269 samples have total relative abundance == 0
    
    # Remove 269 samples from 3467 total Metaphlan profiles
    Metaphlan_profiles_filt <- Metaphlan_profiles[,!SampleFilter_colSums == 0]
    # 3198 samples remaining
    
# Remove samples which failed rel-ab filter from the metadata table
    # find common sample count
    intersectSamples <- intersect(Metadata_Df$SampleID, 
                         colnames(Metaphlan_profiles_filt))
    
    # Filter M4 passed metadata to exclude samples with low abundance
    Metadata_Df <- Metadata_Df %>%
        filter(SampleID %in% intersectSamples)
    # save this table
    write.table(Metadata_Df_sub, 
        "../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv",
        sep = "\t")
    # 3141 samples remaining
    
    # 1. Subset metaphlan sample profiles using new metadata table
    Metaphlan_profiles <- Metaphlan_profiles_filt %>%
        select(all_of(Metadata_Df$SampleID))
        # 3141 samples remaining

    # 2. Subset species which pass the min. abundance criteria
    Species_Total_abundance <- rowSums(Metaphlan_profiles)
    
        # Find number of species with at least 0.001 total abundance
        hist(Species_Total_abundance/3141)    

        table(Species_Total_abundance/3141 > 0.00001, useNA = 'ifany')
            # 812 species
        
        table(Species_Total_abundance/3141 > 0.0001, useNA = 'ifany')
            # 484 species
        
        table(Species_Total_abundance/3141 > 0.001, useNA = 'ifany')
            # 161 species
    
        table(Species_Total_abundance/3141 > 0.01, useNA = 'ifany')
            # 22 species have 10% relative abundance
        
    # Filter species with <0.001 abundance per 10 samples
    # Filter that combines prevalence and abundance cutoffs
    min_samples = 10
    Metaphlan_filtered <- Metaphlan_profiles[rowSums(Metaphlan_profiles >= 0.0001) >= min_samples, ]
        # keeps only species that have a relative abundance of at least 0.0001 in at least 10 samples (min_samples)
        # 1497 species present in at least 10 samples
            # 714 species present with at least 0.0001 abundance in any 10 samples
    
    # sanity check
    Sample_Total_aB <- colSums(Metaphlan_filtered)
        # how many have abundance = 0
    table(Sample_Total_aB > 1)
        View(data.frame(Sample_Total_aB))
        # Some samples have abundance == 0
        
        # which samples have total abundance > 1
        ToRm <- names(which(Sample_Total_aB == 0))
        ToRm <- c("ERR4791684", "SRR12149834", "SRR12397846", "SRR13336638", "S_65306_ID1896_084_S295", "S_69900_ID1896_375_S109")
        '
        [1] "ERR4791684"              "SRR12149834"            
        [3] "SRR12397846"             "SRR13336638"            
        [5] "S_65306_ID1896_084_S295" "S_69900_ID1896_375_S109"
        '
        # check the distribution
        hist(Sample_Total_aB, breaks = 100, main = "Sample total relative abundances",
             xlab = "Total abundance", col = "skyblue")
        abline(v = 0.001, col = "red", lty = 2)
        
    # save this table without the empty samples
    Metaphlan_filtered <- Metaphlan_filtered %>% 
        select(all_of(names(which(Sample_Total_aB > 0))))
        # 3145 samples remaining
    
    write.table(Metaphlan_filtered, 
                "./MiPORT_All_RT_featureTable_sp_filt_1e-4_Relab_per_Min10_batchCorrected_samples.tsv",
                sep = "\t")
    # 3145 samples & 714 species
    
    Metadata_Df_sub <- Metadata_Df %>%
        filter(!SampleID %in% ToRm)

    # save this table
    write.table(Metadata_Df_sub, 
                "../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv",
                sep = "\t")
    