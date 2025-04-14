# Script to filter out samples from global Metaphlan species table according to the RT category

# outputs 3 feature tables with samples from each RT category

# Code to filter Metadata and also filter beta diversity matrix tables
# load libraries
library(vegan)
library(tidyverse)

# read Data
FeatureTable_Df <- read.table("../MetaPhlan4_results/Metaphlan4_species_filtered_Min2_samples.tsv", header = T)

# Get Row and colnames from M4 table
SampleNames <- colnames(FeatureTable_Df[,-1])
RowNames <- row.names(FeatureTable_Df)

# sanity check
print(paste("You have ", ncol(FeatureTable_Df)-1, " samples", nrow(FeatureTable_Df), " species in the feature table"))
# [1] "You have  2935  samples 1497  species in the feature table"

# read Metadata table
Metadata_df <- read.table("../MiPORT_Metadata_downstream_filtered_M4_passed_v2.tsv", header = T, sep = "\t")

print(paste("You have ", nrow(Metadata_df), " number of samples in metadata"))
# [1] "You have  2935  number of samples in metadata"

glimpse(Metadata_df)

# check how many samples are in both tables
RetainSample <- (SampleNames %in% Metadata_df$SampleID)

table(RetainSample)
# 2935 samples present in both

# Filter samples from each RT category and write tables

# overview
table(Metadata_df$RT_category)
    #   IRT  LRT  URT 
    #   720  585 1630 

    # 1. filter LRT
    SamplesToKeep <- Metadata_df %>% 
        filter(RT_category == "LRT") %>%
        pull(SampleID)
    
    # sanity check
    length(SamplesToKeep)
    
    # filter main table
    tempDf <- FeatureTable_Df %>% 
        select(all_of(c("taxonomy", SamplesToKeep)))

    # write this df into file
    write.table(tempDf, "../MetaPhlan4_results/MiPORT_LRT_featureTable_species_filtered_Min2_samples.tsv", sep = '\t', quote = F)

    # 2. filter IRT
    SamplesToKeep <- Metadata_df %>% 
        filter(RT_category == "IRT") %>%
        pull(SampleID)
    
    # sanity check
    length(SamplesToKeep)
    
    # filter main table
    tempDf <- FeatureTable_Df %>% 
        select(all_of(c("taxonomy", SamplesToKeep)))
    
    # write this df into file
    write.table(tempDf, "../MetaPhlan4_results/MiPORT_IRT_featureTable_species_filtered_Min2_samples.tsv", sep = '\t', quote = F)
    
    # 3. filter URT
    SamplesToKeep <- Metadata_df %>% 
        filter(RT_category == "URT") %>%
        pull(SampleID)
    
    # sanity check
    length(SamplesToKeep)
    
    # filter main table
    tempDf <- FeatureTable_Df %>% 
        select(all_of(c("taxonomy", SamplesToKeep)))
    
    # write this df into file
    write.table(tempDf, "../MetaPhlan4_results/MiPORT_URT_featureTable_species_filtered_Min2_samples.tsv", sep = '\t', quote = F)
