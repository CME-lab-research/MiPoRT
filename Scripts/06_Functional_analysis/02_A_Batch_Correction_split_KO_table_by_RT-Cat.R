# Script to split the global feature table according to the RT category
# Input: Base filtered, relative abundance KO (extended) feature table:  02_Humann4_merged_KO_extended_relab.tsv
# Outputs 3 feature tables with samples from each RT category.

# load libraries
library(vegan)
library(tidyverse)

# read Data
FeatureTable_Df <- read.table("../02_Humann4_merged_KO_extended_relab.tsv",
                              header = T,
                              sep = "\t")

# Get Row and colnames from humann KO-table
SampleNames <- colnames(FeatureTable_Df[,-1])
RowNames <- FeatureTable_Df %>% 
    pull(KO_ext_relab)

# logging and sanity check
writeLines(paste("Sample count:", dim(FeatureTable_Df)[2] - 1, 
                 "\nFeature count:", dim(FeatureTable_Df)[1]))
    # Sample count: 3135 
    # Feature count: 9048

# read Metadata table
Metadata_df <- read.table("../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", 
                          header = T, sep = "\t")

print(paste("You have ", nrow(Metadata_df), " number of samples in metadata"))
# [1] "You have  2935  number of samples in metadata"

glimpse(Metadata_df)

# check how many samples are in both tables
RetainSample <- (SampleNames %in% Metadata_df$SampleID)

table(RetainSample)
# 3135 samples present in both

# Filter samples from each RT category and write tables

# overview
table(Metadata_df$RT_category)
#   IRT  LRT  URT 
#   719  578 1838 

# 1. filter LRT
    SamplesToKeep <- Metadata_df %>% 
        filter(RT_category == "LRT") %>%
        pull(SampleID)

    # sanity check
    length(SamplesToKeep) # 578

    # filter main table
    tempDf <- FeatureTable_Df %>% 
        select(all_of(
            c("KO_ext_relab", SamplesToKeep)
            ))

    # write this df into file
    write.table(tempDf, 
                "./MiPORT_LRT_featureTable_KO_ext.tsv", 
                sep = '\t', 
                quote = F)

# 2. filter IRT
    SamplesToKeep <- Metadata_df %>% 
        filter(RT_category == "IRT") %>%
        pull(SampleID)

    # sanity check
    length(SamplesToKeep) # 719

    # filter main table
    tempDf <- FeatureTable_Df %>% 
        select(all_of(c("KO_ext_relab", SamplesToKeep)))
    
    # write this df into file
    write.table(tempDf, 
                "./MiPORT_IRT_featureTable_KO_ext.tsv", 
                sep = '\t', 
                quote = F)

# 3. filter URT
    SamplesToKeep <- Metadata_df %>% 
        filter(RT_category == "URT") %>%
        pull(SampleID)
    
    # sanity check
    length(SamplesToKeep) # 1838
    
    # filter main table
    tempDf <- FeatureTable_Df %>% 
        select(all_of(c("KO_ext_relab", SamplesToKeep)))
    
    # write this df into file
    write.table(tempDf, 
                "./MiPORT_URT_featureTable_KO_ext.tsv", 
                sep = '\t', 
                quote = F)

        
    
    
    
    
    
    
