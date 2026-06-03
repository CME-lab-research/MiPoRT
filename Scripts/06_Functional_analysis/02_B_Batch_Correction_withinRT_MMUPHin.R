library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)

# Read metadata file
  Metadata_df <- read.table(
    "../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv",
    header = T, 
    sep = "\t", 
    check.names = F)

  # remove blank rows
  Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]
  Metadata_df <- Metadata_df %>% filter(SampleType != 'Anterior_nares')

# MMUPHin fn data ingetrity check
check_feature_abd <- function(feature_abd) {
  #print(colSums(feature_abd))
  # Errors if has missing values
  if(any(is.na(feature_abd)))
    stop("Found missing values in the feature table!")
  # Errors if has negative values
  # if statement condition is okay because NA's would've been checked
  if(any(feature_abd < 0))
    stop("Found negative values in the feature table!")
  # Returns "proportions" if all values are less than or equal to one
  # Non-negativeness has been checked by previous if statement
  if(all(feature_abd <= 1)) {
    return("proportions")
  }
  # Returns "counts" if all values are integers
  else if(all(feature_abd == floor(feature_abd))) {
    return("counts")
  }
  # Errors if is neither proportions nor counts (i.e., has values greater than
  # one that are not integers)
  else
    stop("Feature table does not appear to be either proportions or counts!")
}

######################################################################################
##################### Batch correction for IRT Sample Type ###########################
######################################################################################

# Read Humann-4 profile file for IRT category
Humann_profiles <- read.table("./MiPORT_IRT_featureTable_KO_ext.tsv", 
                              header = T, 
                              sep = "\t", 
                              fill = T, 
                              check.names = F)

# logging and sanity check
writeLines(paste("Sample count:", dim(Humann_profiles)[2] - 1, 
                 "\nFeature count:", dim(Humann_profiles)[1]))
  # Sample count:   719 
  # Feature count: 9048

  # Save sample IDs
  Samples_Total_before <- colSums(Humann_profiles[,-1])
  table(Samples_Total_before)

  # Following samples have 0 abundances; Probably only had singletons
  SamplestoRm <- names(which(Samples_Total_before == 0))
  names(which(Samples_Total_before== 0))
  
  # If empty, good then move on.
  # If not the rm failed samples.
  Humann_profiles <- Humann_profiles %>%
    select(c("taxonomy",!all_of((SamplestoRm)))) # breaks code
  
  # subset Metadata df
  Metadata_df_sub <- Metadata_df %>%
    filter(RT_category == "IRT")

# Step 1: Data preparation (Common for all sample types)
  # Get a list of Humann-4 profile samples to retain
  retainSamples <- intersect(
    (colnames(Humann_profiles)[-1]), 
    Metadata_df_sub$SampleID)
  
  # Check how many are common between both tables
  length(retainSamples)
  #RM_retainSamples <- intersect((colnames(Humann_profiles)[-1]), Metadata_df_RM$SampleID)
  
  # Prepare feature abundances: To be provided as a feature-by-sample matrix
  row.names(Humann_profiles) <- Humann_profiles$KO_ext_relab
  Humann_profiles <- Humann_profiles[retainSamples]

  # logging and sanity check
  writeLines(paste("Sample count:", dim(Humann_profiles)[2], 
                   "\nFeature count:", dim(Humann_profiles)[1]))
  # Sample count:   719 
  # Feature count: 9048

  # subset metadata df to match
  Metadata_df_sub <- Metadata_df_sub %>% 
    filter(SampleID %in% retainSamples)
  row.names(Metadata_df_sub) <- Metadata_df_sub$SampleID

  # sanity check
  dim(Metadata_df_sub)

# Step 2: Numeric values formatting
  # Convert all columns of a data frame to numeric or integer
  Humann_profiles[] <- lapply(Humann_profiles, function(x) {
    if (is.factor(x)) {
      as.numeric(as.character(x))  # Convert factors to numeric
    } else if (is.character(x)) {
      as.numeric(x)  # Convert character to numeric
    } else {
      as.numeric(x)  # Ensure numeric columns remain numeric
    }
  })

# Step 3: data table checks from MMUPHIN (Common for all sample types)

  # Check if cols are numeric values
  Sum_abundances <- colSums(Humann_profiles)
  check_feature_abd(Humann_profiles) 
  # "proportions" 
  # Passed :)

  # Sanity check batches
  table(Metadata_df_sub$BioProject, Metadata_df_sub$SampleType)
  table(Metadata_df_sub$SampleType)

# Check if we match sample names between metadata and profiles
# Step 3.a: Ensure Metadata_df_sub row names match colnames of Humann_profiles_scaled
  Metadata_df_sub <- Metadata_df_sub[
    match(
      colnames(Humann_profiles), 
      row.names(Metadata_df_sub)
      ), 
    ]

# Step 3.b: Check if order is now identical
  identical(colnames(Humann_profiles), 
            row.names(Metadata_df_sub))  
  # Should return TRUE

# add factors
  Metadata_df_sub$BioProject <- 
    factor(Metadata_df_sub$BioProject)

  # sanity check
  table(Metadata_df_sub$BioProject)
  table(Metadata_df_sub$SampleType)
  levels(Metadata_df_sub$BioProject)

# Make a copy of the df
Humann_profiles_cp <- Humann_profiles[,Metadata_df_sub$SampleID] 

# Step 4: Batch correction
# For IRT we subset for Sputum samples and correct only for that since no other samples types are present in multiple Bio-Projects
# Select sampleType
  sample = "Sputum"

  # Subset again metadata and feature table for the current sample type
  subset_metadata <- Metadata_df_sub[
    Metadata_df_sub$SampleType == sample, ]
  subset_features <- Humann_profiles[, subset_metadata$SampleID]

  # Drop levels
  subset_metadata$BioProject <- droplevels(subset_metadata$BioProject)

# Run MMUPHin batch correction for this sample type
# Adjust for batch effects while controlling for the effect of SampleType
  Simple_fit_Humann_profiles <- adjust_batch(
    feature_abd = subset_features,  # Taxa proportions
    batch = "BioProject",        # Batch effect to correct
    data = subset_metadata, # Metadata file
    control = list(verbose = FALSE)
  )
  
  # Get the corrected abundance table now
  Humann_profiles_scaled_adj <- 
    data.frame(Simple_fit_Humann_profiles$feature_abd_adj)

  # Merge new Sputum sample abundances with other sampleType abundances
  # get sampleIds which are not Sputum to join below
  SamplesIDs_toRetain <- Metadata_df_sub %>% 
    filter(SampleType != "Sputum") %>%
    pull(SampleID)

  # Join Sputum and other IRT sample types
  Humann_profiles_scaled_adj <- cbind(
    Humann_profiles[, SamplesIDs_toRetain], 
    data.frame(Simple_fit_Humann_profiles$feature_abd_adj))

# Step 5: Save results
  write.table(Humann_profiles_scaled_adj, 
              "MiPORT_BC-IRT_featureTable_KO_ext.tsv", 
              row.names = T, 
              sep = "\t", 
              quote = F)

# Step 6: Check effect reduction
# Test variance explained now
  library(vegan, quietly = TRUE)
  set.seed(8284)

  # random check
  Samples_Total_before <- colSums(Humann_profiles)
  Samples_Total_after <- colSums(Humann_profiles_scaled_adj)

  # calc BC dist
  D_before <- vegdist(t(Humann_profiles))
  D_after <- vegdist(t(Humann_profiles_scaled_adj))

  # Calculate R-sq with Adonis fit to predict Bio-project
  fit_adonis_before <- adonis2(D_before ~ BioProject, 
                               data = Metadata_df_sub,
                               strata = Metadata_df_sub$SampleType, 
                               parallel = 8) 
  
  fit_adonis_after <- adonis2(D_after ~ BioProject, 
                              data = Metadata_df_sub,
                              strata = Metadata_df_sub$SampleType,  
                              parallel = 8) 

  # Check effect size reduction
  print(fit_adonis_before)
  print(fit_adonis_after)

  # Output summarised in a table below.
  # MMUPHin reduced BioProject batch effects in KO profiles from ~30% to ~22% variance explained, indicating effective but partial correction.
  "
  Metric	Before correction	After correction	Change
  Df (BioProject)	9 9	—
  R² (BioProject effect)	0.298	0.218	↓ −0.080 (~27% reduction)
  F-statistic	33.45	21.94	↓ −34%
  p-value	0.001 ***	0.001 ***	unchanged (still significant)
  Residual R²	0.702	0.782	↑ +0.080
  Total variance	1.000	1.000	—
  "
######################################################################################
##################### Batch correction for LRT Sample Type ###########################
######################################################################################

# Read filtered Humann-4 profile file for LRT category
  Humann_profiles <- 
    read.table("MiPORT_LRT_featureTable_KO_ext.tsv", 
                              header = T, 
                              sep = "\t", 
                              fill = T, 
                              check.names = F)

# Save sample IDs
  Samples_Total_before <- colSums(Humann_profiles[,-1])
  table(Samples_Total_before)

# Rm samples with Total abundance 0
# The following samples have 0 abundances; Probably only had singletons
  SamplestoRm <- names(which(Samples_Total_before == 0))
  names(which(Samples_Total_before== 0))

# rm failed samples: Mv this up later
# Humann_profiles <- Humann_profiles %>% select(c("taxonomy",!all_of((SamplestoRm)))) # breaks code?
# if none found then proceed without filtering Samples out

# subset Metadata df
  Metadata_df_sub <- Metadata_df %>%
    filter(RT_category == "LRT")

# Step 1: Data preparation (Common for all sample types)
  # Get a list of Humann-4 profile samples to retain
  retainSamples <- intersect(
    (colnames(Humann_profiles)[-1]), 
    Metadata_df_sub$SampleID
    )
  
  # Check how many are common between both tables
  length(retainSamples)

# Prepare feature abundances be provided as a feature-by-sample matrix
  row.names(Humann_profiles) <- Humann_profiles$KO_ext_relab
  Humann_profiles <- Humann_profiles[retainSamples]
  
  # sanity check
  dim(Humann_profiles)

  # subset metadata df to match
    Metadata_df_sub <- Metadata_df_sub %>% 
      filter(SampleID %in% retainSamples)
    
    row.names(Metadata_df_sub) <- Metadata_df_sub$SampleID
  
    # sanity check
    dim(Metadata_df_sub)

# Step 2: Numeric values formatting
  # Convert all columns of a data frame to numeric or integer
  Humann_profiles[] <- lapply(Humann_profiles, function(x) {
    if (is.factor(x)) {
      as.numeric(as.character(x))  # Convert factors to numeric
    } else if (is.character(x)) {
      as.numeric(x)  # Convert character to numeric
    } else {
      as.numeric(x)  # Ensure numeric columns remain numeric
    }
  })

# Step 3: data table checks from MMUPHIN (Common for all sample types)

  # Check if cols are numeric values
  Sum_abundances <- colSums(Humann_profiles)
  check_feature_abd(Humann_profiles) 
  # "proportions" 
  # Passed :)

  # Sanity check batches
  table(Metadata_df_sub$BioProject, Metadata_df_sub$SampleType)
  table(Metadata_df_sub$SampleType)

# Check if we match sample names between metadata and profiles
  # Step 3.a: Ensure Metadata_df_sub row names match colnames of Humann_profiles_scaled
  Metadata_df_sub <- 
    Metadata_df_sub[
      match(colnames(Humann_profiles), 
            row.names(Metadata_df_sub)), 
      ]

  # Step 3.b: Check if order is now identical
  identical(colnames(Humann_profiles), 
            row.names(Metadata_df_sub))  
  # Should return TRUE

  # add factors
  Metadata_df_sub$BioProject <- 
    factor(Metadata_df_sub$BioProject)

  # sanity check
  table(Metadata_df_sub$BioProject)
  table(Metadata_df_sub$SampleType)
  levels(Metadata_df_sub$BioProject)

# Step 4: Batch correction
  # For single sampleType we don't subset 

  # Run MMUPHin batch correction for this sample type
  Simple_fit_Humann_profiles <- 
    adjust_batch(
    feature_abd = Humann_profiles,  # KO-ext. proportions
    batch = "BioProject",        # Batch effect to correct
    data = Metadata_df_sub,     # Metadata file
    control = list(verbose = FALSE)
  )
  
  # Get the corrected abundance table now
  Humann_profiles_scaled_adj <- data.frame(Simple_fit_Humann_profiles$feature_abd_adj)

  # Step 5: Save results
  write.table(Humann_profiles_scaled_adj, 
              "MiPORT_BC-LRT_featureTable_KO_ext.tsv", 
              row.names = T, 
              sep = "\t", 
              quote = F)

# Step 6: Check effect reduction
  # Test variance explained now
  #library(vegan, quietly = TRUE)
  #set.seed(8284)

  # random check
  Samples_Total_before <- colSums(Humann_profiles)
  Samples_Total_after <- colSums(Humann_profiles_scaled_adj)

  # calc BC dist
  D_before <- vegdist(t(Humann_profiles))
  D_after <- vegdist(t(Humann_profiles_scaled_adj))

# Calculate R-sq with Adonis fit to predict Bio-project
  fit_adonis_before <- adonis2(D_before ~ BioProject, 
                               data = Metadata_df_sub, 
                               parallel = 8) 
  # For LRT remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

  fit_adonis_after <- adonis2(D_after ~ BioProject, 
                              data = Metadata_df_sub, 
                              parallel = 8) 
  # For LRT remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

# Check effect size reduction
print(fit_adonis_before)
print(fit_adonis_after)


######################################################################################
##################### Batch correction for URT Sample Type ###########################
######################################################################################


# Read Humann-4 profile file for URT category
Humann_profiles <- read.table("./MiPORT_URT_featureTable_KO_ext.tsv", 
                              header = T, 
                              sep = "\t", 
                              fill = T, 
                              check.names = F)

# NOTE: subset for Buccal_mucosa & Saliva samples since these are the only 2 with multiple batches
  
  # sanity check
  Samples_Total_before <- round(colSums(Humann_profiles[,-1]), 1)
  table(Samples_Total_before)

  # Following samples have 0 abundances; Probably only had singletons
  SamplestoRm <- names(which(Samples_Total_before == 0))
  names(which(Samples_Total_before== 0))

  # subset Metadata df
  Metadata_df_sub <- Metadata_df %>%
    filter(RT_category == "URT")

# Step 1: Data preparation (common for all sample types)
  
  # Get a list of Humann-4 profile samples to retain
  retainSamples <- intersect((colnames(Humann_profiles)[-1]), Metadata_df_sub$SampleID)
  # Check how many are common between both tables
  length(retainSamples)

# Prepare feature abundances be provided as a feature-by-sample matrix
  row.names(Humann_profiles) <- Humann_profiles$KO_ext_relab
  Humann_profiles <- Humann_profiles[retainSamples]

  # sanity check
  dim(Humann_profiles)

  # subset metadata df to match
  Metadata_df_sub <- Metadata_df_sub %>% 
    filter(SampleID %in% retainSamples)
  
  row.names(Metadata_df_sub) <- 
    Metadata_df_sub$SampleID

  # sanity check
  dim(Metadata_df_sub)

# Step 2: Numeric values formatting
  # Convert all columns of a data frame to numeric or integer
  Humann_profiles[] <- lapply(Humann_profiles, function(x) {
    if (is.factor(x)) {
      as.numeric(as.character(x))  # Convert factors to numeric
    } else if (is.character(x)) {
      as.numeric(x)  # Convert character to numeric
    } else {
      as.numeric(x)  # Ensure numeric columns remain numeric
    }
  })

# Step 3: data table checks from MMUPHIN (Common for all sample types)
  
  # Check if cols are numeric values
  Sum_abundances <- colSums(Humann_profiles)
  
  # MMUPHN check
  check_feature_abd(Humann_profiles) 
  # "proportions"
  # PASSED check. Proceed.

  # Sanity check batches
  table(Metadata_df_sub$BioProject, 
        Metadata_df_sub$SampleType)
  
  table(Metadata_df_sub$SampleType)
  
  # Subset for sample type with multiple batches
  # Keep only these 2 sample Types as they have different batches
  SampleTypes_to_Correct <- c("Buccal_mucosa", "Saliva")
  Metadata_df_sub_URT <- Metadata_df_sub %>%
    filter(SampleType %in% c(SampleTypes_to_Correct))

  # sanity check
  table(Metadata_df_sub_URT$SampleType)

  # Get a list of Humann-4 profile samples to retain
  retainSamples <- intersect(
    (colnames(Humann_profiles)), 
    Metadata_df_sub_URT$SampleID)
  
  length(retainSamples)
  dim(Humann_profiles) # before 9048 1630

  # subset Humann-4 profiles
  Humann_profiles_scaled_sub <- Humann_profiles[retainSamples]
  dim(Humann_profiles_scaled_sub) # after 9048 1342 (1027 + 315)

# Check if we match sample names between metadata and profiles
  # Step 3.a: Ensure Metadata_df_sub_URT row names match colnames of Humann_profiles_scaled_sub
  Metadata_df_sub_URT <- 
    Metadata_df_sub_URT[
      match(colnames(Humann_profiles_scaled_sub), 
            row.names(Metadata_df_sub_URT)), 
      ]

  # Step 3.b: Check if order is now identical
  identical(colnames(Humann_profiles_scaled_sub), row.names(Metadata_df_sub_URT))  
  # Should return TRUE

  # add factors
  Metadata_df_sub_URT$BioProject <- 
    factor(Metadata_df_sub_URT$BioProject)

  # sanity check
  table(Metadata_df_sub_URT$BioProject)
  table(Metadata_df_sub_URT$SampleType)

# Step 4: Batch correction
  # URT has 2 sampletypes to correct. We subset and then correct for each

  # Make a copy of the input URT df
  Humann_profiles_scaled_sub_OG <- Humann_profiles_scaled_sub
  Metadata_df_sub_URT_OG <- Metadata_df_sub_URT

# A. Correct for Saliva 
  # subset for 1 sampleType in this
  SamplesIDs_toRetain <- Metadata_df_sub_URT_OG %>% 
    filter(SampleType == "Saliva") %>% # Either Buccal_mucosa or Saliva
    pull(SampleID)

  # sanity check
  length(SamplesIDs_toRetain)

  # subset metadata and humann table for these
  Humann_profiles_scaled_sub <- 
    Humann_profiles_scaled_sub_OG[, SamplesIDs_toRetain]
  Metadata_df_sub_URT <- Metadata_df_sub_URT_OG %>% 
    filter(SampleType == "Saliva") 

  # drop unused factors
  Metadata_df_sub_URT$BioProject <- droplevels(Metadata_df_sub_URT$BioProject)

  # sanity check
  dim(Humann_profiles_scaled_sub)
  dim(Metadata_df_sub_URT) 
  table(Metadata_df_sub_URT$BioProject)

  # Run MMUPHin batch correction for this sample type
  Simple_fit_Humann_profiles <- adjust_batch(
    feature_abd = Humann_profiles_scaled_sub,  # Humann proportions
    batch = "BioProject",        # Batch effect to correct
    data = Metadata_df_sub_URT, # Metadata file
    control = list(verbose = FALSE)
  )
  
  # Get the corrected abundance table for saliva
  Humann_profiles_scaled_adj_Saliva <- 
    data.frame(Simple_fit_Humann_profiles$feature_abd_adj)

# B. Correct for Buccal_mucosa 
  # subset for 1 sampleType in this
  SamplesIDs_toRetain <- Metadata_df_sub_URT_OG %>% 
    filter(SampleType == "Buccal_mucosa") %>% # Either Buccal_mucosa or Saliva
    pull(SampleID)

  # sanity check
  length(SamplesIDs_toRetain)

  # subset metadata and metaphlan table for these
  Humann_profiles_scaled_sub <- 
    Humann_profiles_scaled_sub_OG[,SamplesIDs_toRetain]
  
  Metadata_df_sub_URT <- Metadata_df_sub_URT_OG %>% 
    filter(SampleType == "Buccal_mucosa") 

  # drop unused factors
  Metadata_df_sub_URT$BioProject <- droplevels(Metadata_df_sub_URT$BioProject)

  # sanity check
  dim(Humann_profiles_scaled_sub)
  dim(Metadata_df_sub_URT) 
  table(Metadata_df_sub_URT$BioProject)

# Run MMUPHin batch correction for this sample type
Simple_fit_Humann_profiles <- adjust_batch(
  feature_abd = Humann_profiles_scaled_sub,  # humann proportions
  batch = "BioProject",        # Batch effect to correct
  data = Metadata_df_sub_URT, # Metadata file
  control = list(verbose = FALSE)
)

# Get the corrected abundance table for Buccal mucosa
Humann_profiles_scaled_adj_Buccal_mucosa <- 
  data.frame(Simple_fit_Humann_profiles$feature_abd_adj)


# c. Finally, we merge all 3 tables together 
# Humann_profiles_scaled_adj_Buccal_mucosa,
# Humann_profiles_scaled_adj_Saliva, 
# OG_abundance from other sampleTypes (other than Buccal and saliva)

# quick check
table(Metadata_df_sub$SampleType)

# get sampleIds which are not in SampleTypes_to_Correct
SamplesIDs_toRetain <- Metadata_df_sub %>% 
  filter(!(SampleType %in% SampleTypes_to_Correct)) %>%
  pull(SampleID)

length(SamplesIDs_toRetain) # Should be 288 (27 + 15 + 246)

  # Merge
  Humann_profiles_scaled_adj <- 
    cbind(
      Humann_profiles[, SamplesIDs_toRetain],
      Humann_profiles_scaled_adj_Buccal_mucosa,
      Humann_profiles_scaled_adj_Saliva
      )
  
  dim(Humann_profiles_scaled_adj)
  dim(Humann_profiles)

# Step 5: Save results
write.table(Humann_profiles_scaled_adj, 
            "MiPORT_BC-URT_featureTable_KO_ext.tsv", 
            row.names = T, 
            sep = "\t", 
            quote = F)

# Step 6: Check effect reduction

# Test variance explained now
  
  # A. Buccal mucosa adjustments
  dim(Humann_profiles_scaled_sub)
  dim(Humann_profiles_scaled_adj_Buccal_mucosa)
  
  Samples_Total_before <- round(
    colSums(Humann_profiles_scaled_sub), 1)
  Samples_Total_after <- round(
    colSums(Humann_profiles_scaled_adj_Buccal_mucosa), 1)

  table(Samples_Total_after)
  table(Samples_Total_before)

  # calc BC dist
  D_before <- vegdist(t(Humann_profiles_scaled_sub))
  D_after <- vegdist(t(Humann_profiles_scaled_adj_Buccal_mucosa))

  # Calculate R2 with Adonis fit to predict Bio-project
  fit_adonis_before <- adonis2(D_before ~ BioProject, 
                               data = Metadata_df_sub_URT, 
                               parallel = 8) 
    # Add 'strata = Metadata_df_sub$SampleType' param 

  fit_adonis_after <- adonis2(D_after ~ BioProject, 
                              data = Metadata_df_sub_URT, 
                              parallel = 8) 

  print(fit_adonis_before)
  print(fit_adonis_after)

  # R sq reduced from 0.14 to 0.0046
  
  # B. Saliva adjustments
  dim(Humann_profiles_scaled_adj_Saliva)
  dim(Metadata_df_sub_URT)

  Samples_Total_before <- round(
    colSums(Humann_profiles_scaled_sub), 1)
  Samples_Total_after <- 
    round(
      colSums(Humann_profiles_scaled_adj_Saliva), 1)

  table(Samples_Total_before)
  table(Samples_Total_after)

  
  D_before <- vegdist(t(Humann_profiles_scaled_sub))
  D_after <- vegdist(t(Humann_profiles_scaled_adj_Saliva))

# Calculate R2 with Adonis fit to predict Bio-project
  fit_adonis_before <- 
    adonis2(D_before ~ BioProject, 
            data = Metadata_df_sub_URT, 
            parallel = 8) 
  
  fit_adonis_after <- 
    adonis2(D_after ~ BioProject, 
            data = Metadata_df_sub_URT, 
            parallel = 8) 

  print(fit_adonis_before)
  print(fit_adonis_after)
  
  # Saliva: R sq. reduced from 0.25 to 0.019
  
# C. Global URT
  dim(Humann_profiles)
  dim(Humann_profiles_scaled_adj)
  
  Samples_Total_before <- round(colSums(Humann_profiles), 1)
  Samples_Total_after <- round(colSums(Humann_profiles_scaled_adj), 1)
  
  table(Samples_Total_after)
  table(Samples_Total_before)

  # calc BC dist
  D_before <- vegdist(t(Humann_profiles))
  D_after <- vegdist(t(Humann_profiles_scaled_adj))

# Calculate R2 with Adonis fit to predict Bio-project
  # Add 'strata = Metadata_df_sub$SampleType' param 
  fit_adonis_before <- 
    adonis2(D_before ~ BioProject, 
            data = Metadata_df_sub, 
            parallel = 8, 
            strata = Metadata_df_sub$SampleType) 

  fit_adonis_after <- 
    adonis2(D_after ~ BioProject, 
            data = Metadata_df_sub, 
            parallel = 8, 
            strata = Metadata_df_sub$SampleType) 

  print(fit_adonis_before)
  print(fit_adonis_after)

# ---------------------- DONE ---------------------- #

| --------- | --------- | ----- | ----- |
| Metric    | URT       | IRT   | LRT   |
| --------- | --------- | ----- | ----- |
| R² before | 0.272     | 0.298 | 0.211 |
| R² after  | **0.057** | 0.218 | 0.145 |
| Reduction | **~79%**  | ~27%  | ~31%  |
| --------- | --------- | ----- | ----- |

