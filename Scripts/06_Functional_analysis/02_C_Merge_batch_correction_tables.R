# Merge batch corrected samples from the following files
  # LRT: MiPORT_LRT_featureTable_KO_ext.tsv
  # IRT: MiPORT_IRT_featureTable_KO_ext.tsv
  # URT: MiPORT_BC-URT_featureTable_KO_ext.tsv

# load libraries
library(tidyverse)

# read files
  URT_df <- 
    read.table(
      "MiPORT_BC-URT_featureTable_KO_ext.tsv",
      header = T,
      sep = "\t",
      check.names = F)
  
  IRT_df <- 
    read.table("./MiPORT_IRT_featureTable_KO_ext.tsv",
               header = T,
               sep = "\t",
               check.names = F)
  
  # mv KO_ext_relab column to rownames
  row.names(IRT_df) <- IRT_df$KO_ext_relab
  IRT_df <- IRT_df[,-1]
  
  # Read LRT table
  LRT_df <- 
    read.table("./MiPORT_LRT_featureTable_KO_ext.tsv",
               header = T,
               sep = "\t",
               check.names = F)
  
  # mv KO_ext_relab column to rownames
  row.names(LRT_df) <- LRT_df$KO_ext_relab
  LRT_df <- LRT_df[,-1]
  
  # sanity check
  hist(colSums(URT_df))
  hist(colSums(IRT_df))
  hist(colSums(LRT_df))
  
  # Anterior nares
  Other_samples_df <- read.table(
    "../02_Humann4_merged_KO_extended_relab.tsv",
    header = T,
    sep = '\t',
    check.names = F
  )
  

  hist(colSums(Other_samples_df[,-1]))

# URT samples were batch corrected individually for sample types Buccal_mucosa & saliva
# IRT samples from Sputum were batch corrected
# LRT samples from BAL were batch corrected

# column bind 3 dfs (URT, IRT, LRT df)
  Merged_df <- cbind(URT_df, 
                     IRT_df, 
                     LRT_df)
  # 720 + 579 + 1631 = 2930
  
  Merged_df <- Merged_df %>%
    rownames_to_column(var="KO_ext_relab")
  
# save batch corrected global table
  write.table(Merged_df,
              "./03_BC_Humann4_merged_KO_extended_relab.tsv",
              sep = '\t',
              row.names = F,
              quote = F)
  
  # The above file was without Anterior nares. 
  
  # Save a file with Anterior nare samples also
    '
                  IRT  LRT  URT
  Anterior_nares    0    0  208
  BAL               0  578    0
  Buccal_mucosa     0    0 1027
  Nasal_Swab        0    0   27
  Oral_swab         0    0  246
  Other            49    0   15
  Saliva            0    0  315
  Sputum          193    0    0
  Supraglottal     59    0    0
  Tongue_dorsum   418    0    0
    '
    library(dplyr)
    
    # Find sample columns present in Other_samples_df but missing in Merged_df
    Distinct_samples <- setdiff(colnames(Other_samples_df), colnames(Merged_df))
    
    # Remove feature ID column if it gets picked up
    Distinct_samples <- setdiff(Distinct_samples, "KO_ext_relab")
    
    length(Distinct_samples)
    # 208 samples are missing from Merged_df
    
    # Add missing sample columns by matching on KO_ext_relab
    Merged_final <- Merged_df %>%
      left_join(
        Other_samples_df %>%
          select(KO_ext_relab, all_of(Distinct_samples)),
        by = "KO_ext_relab"
      )
    
    # Check all required columns exist
    if (!all(colnames(Other_samples_df) %in% colnames(Merged_final))) {
      missing_cols <- setdiff(colnames(Other_samples_df), colnames(Merged_final))
      stop(paste("Missing columns in Merged_final:", paste(missing_cols, collapse = ", ")))
    }
    
    # Reorder columns to match Other_samples_df exactly
    Merged_final <- Merged_final[, colnames(Other_samples_df), drop = FALSE]
    
    # Final sanity check
    identical(colnames(Merged_final), colnames(Other_samples_df))
    # TRUE
    
    # Other sanity checks
    anyDuplicated(Merged_df$KO_ext_relab)
    anyDuplicated(Other_samples_df$KO_ext_relab)
 
    # Save this final table for merged batch corrected filtered profiles
    # from metaphlan
    write.table(Merged_final,
                "./MiPORT_BC-RT_withAN_featureTable_KO_ext.tsv",
                sep = '\t',
                row.names = F,
                quote = F)
  
  
