'Script to calculate diversity metrics in R post Batch correction.
Input: Metaphlan species abundance feature table (Batch corrected values)
    MiPORT_RT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv

'

# ---- Load Required Libraries ----
library(tidyverse)
library(optparse)
source("./calculate_diversity.R")

# Read M4 table
Metaphlan_BC_profiles <- file.path(
  "../04_MMuphin_Batch_correction/MiPORT_All_RT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv")

if(file.exists(Metaphlan_BC_profiles)){
  cat("Proceeding with diversity calculations")
}else{(
  cat("File not found")
)}

# Function to calculate a set of alpha diversities
calc_alpha  <- function(feature_table) {
  cat("Function to calculate alpha diversity from Metaphlan profiles.\n")
  cat("Date:", format(Sys.time(), "%d-%b-%H_%M"), "\n\n")
  
  # ---- Define Input Paths ----
  ProjectID <- 'MiPoRT_all_RT'
  Date_Time <- format(Sys.time(), "%d-%b-%H_%M")
  WorkDir <- "C:/Users/o_shinde/Nextcloud/AG Moissl-Eichinger_people/Tejus/MiPoRT Paper/MiPORT_DA/04_calculate_diversity_postBC" # to save output file
  LogDir <- file.path(WorkDir, "Logs")
  Merged_profiles <- Metaphlan_BC_profiles
  
  # ---- Define Diversity Metrics ----
  alphaDiv_options <- c("richness", "shannon", "simpson", "gini")
  
  # Loop through metrics and call external R script
    for (alphaDiv in alphaDiv_options) {
      cmd <- sprintf("Rscript %s/../Scripts/calculate_diversity.R -f %s -d alpha -m %s -o %s",
                     WorkDir, Merged_profiles, alphaDiv, WorkDir)
      cat("Running command:\n", cmd, "\n")
      # calc
      calculate_diversity(profile_file = Merged_profiles, 
                          diversity_type = "alpha", 
                          metric = alphaDiv)
    }
  
  cat("End program\n")
  
  
}
