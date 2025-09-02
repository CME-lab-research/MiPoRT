library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)

# Read metadata file
  Metadata_df <- read.table("../MiPORT_Metadata_downstream_filtered_M4_passed.tsv", header = T, sep = "\t", check.names = F)

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

# Read M4 profile file for IRT category
Metaphlan_profiles <- read.table("../MetaPhlan4_results/RT_based_feature_table_subsets/MiPORT_IRT_featureTable_species_filtered_Min2_samples.tsv", 
  header = T, 
  sep = "\t", 
  fill = T, 
  check.names = F)

  # Save sample IDs
  Samples_Total_before <- colSums(Metaphlan_profiles[,-1])
  table(Samples_Total_before)

  # Following samples have 0 abundances; Probably only had singletons
  SamplestoRm <- names(which(Samples_Total_before == 0))
  names(which(Samples_Total_before== 0))

  # rm failed samples: Mv this up later
  Metaphlan_profiles <- Metaphlan_profiles %>%
    select(c("taxonomy",!all_of((SamplestoRm)))) # breaks code?
  # if none found then proceed without filtering Samples out

  # subset Metadata df
  Metadata_df_sub <- Metadata_df %>%
    filter(RT_category == "IRT")

# Step 1: Data preparation (Common for all sample types)
# Get a list of M4 profile samples to retain
retainSamples <- intersect((colnames(Metaphlan_profiles)[-1]), Metadata_df_sub$SampleID)
  # Check how many are common between both tables
 length(retainSamples)
#RM_retainSamples <- intersect((colnames(Metaphlan_profiles)[-1]), Metadata_df_RM$SampleID)

# Prepare feature abundances be provided as a feature-by-sample matrix
row.names(Metaphlan_profiles) <- Metaphlan_profiles$taxonomy
Metaphlan_profiles <- Metaphlan_profiles[retainSamples]

  # sanity check
  dim(Metaphlan_profiles)

# subset metadata df to match
Metadata_df_sub <- Metadata_df_sub %>% filter(SampleID %in% retainSamples)
row.names(Metadata_df_sub) <- Metadata_df_sub$SampleID

  # sanity check
  dim(Metadata_df_sub)

# Step 2: Numeric values formatting
# Convert all columns of a data frame to numeric or integer
Metaphlan_profiles[] <- lapply(Metaphlan_profiles, function(x) {
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
    Sum_abundances <- colSums(Metaphlan_profiles)
    check_feature_abd(Metaphlan_profiles) 
    # will FAIL because MMUPHIN expects a different scale of sample values. The total sum per sample is 100. Therefore, we scale it to 1.

  # Scaling the data from 0-100 range to 0-1
    Metaphlan_profiles_scaled <- Metaphlan_profiles / 100
    check_feature_abd(Metaphlan_profiles_scaled)
    
  # Sanity check batches
    table(Metadata_df_sub$BioProject, Metadata_df_sub$SampleType)
    table(Metadata_df_sub$SampleType)

  # Check if we match sample names between metadata and profiles
  # Step 3.a: Ensure Metadata_df_sub row names match colnames of Metaphlan_profiles_scaled
  Metadata_df_sub <- Metadata_df_sub[match(colnames(Metaphlan_profiles_scaled), row.names(Metadata_df_sub)), ]

  # Step 3.b: Check if order is now identical
  identical(colnames(Metaphlan_profiles_scaled), row.names(Metadata_df_sub))  
  # Should return TRUE

  # add factors
  Metadata_df_sub$BioProject <- factor(Metadata_df_sub$BioProject)

  # sanity check
  table(Metadata_df_sub$BioProject)
  table(Metadata_df_sub$SampleType)
  levels(Metadata_df_sub$BioProject)

  # Make a copy of the df
  Metaphlan_profiles_scaled_sub <- Metaphlan_profiles_scaled[,Metadata_df_sub$SampleID] 

# Step 4: Batch correction
  # For IRT we subset for Sputum samples and correct only for that since no other samples types are present in multiple BioProjects
  # Select sampleType
  sample = "Sputum"

  # Subset again metadata and feature table for the current sample type
  subset_metadata <- Metadata_df_sub[Metadata_df_sub$SampleType == sample, ]
  subset_features <- Metaphlan_profiles_scaled[, subset_metadata$SampleID]
  
  # Drop levels
  subset_metadata$BioProject <- droplevels(subset_metadata$BioProject)
  
  # Run MMUPHin batch correction for this sample type
  # Adjust for batch effects while controlling for the effect of SampleType
  Simple_fit_Metaphlan_profiles <- adjust_batch(
    feature_abd = subset_features,  # Taxa abundance table scaled to 1
    batch = "BioProject",        # Batch effect to correct
    #covariates = c("Healthy"),
    data = subset_metadata, # Metadata file
    control = list(verbose = FALSE)
  )
  
  # Get the corrected abundance table now
  Metaphlan_profiles_scaled_adj <- data.frame(corrected_data$feature_abd_adj)
  
  # Merge new Sputum sample abundances with other sampleType abundances
  # get sampleIds which are not Sputum to join below
  SamplesIDs_toRetain <- Metadata_df_sub %>% 
    filter(SampleType != "Sputum") %>%
    pull(SampleID)
  
  # Join Sputum and other IRT sample types
  Metaphlan_profiles_scaled_adj <- cbind(Metaphlan_profiles_scaled[, SamplesIDs_toRetain], data.frame(corrected_data$feature_abd_adj))

# Step 5: Save results
  write.table(Metaphlan_profiles_scaled_adj, "../MetaPhlan4_results/MiPORT_IRT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", row.names = T, sep = "\t", quote = F)
  
# Step 6: Check effect reduction
  # Test variance explained now
  library(vegan, quietly = TRUE)
  set.seed(8284)

  # random check
  Samples_Total_before <- colSums(Metaphlan_profiles_scaled)
  Samples_Total_after <- colSums(Metaphlan_profiles_scaled_adj)

  # calc BC dist
  D_before <- vegdist(t(Metaphlan_profiles_scaled))
  D_after <- vegdist(t(Metaphlan_profiles_scaled_adj))

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

######################################################################################
##################### Batch correction for LRT Sample Type ###########################
######################################################################################

  # Read filtered M4 profile file for LRT category
  Metaphlan_profiles <- read.table("../MetaPhlan4_results/RT_based_feature_table_subsets/MiPORT_LRT_featureTable_species_filtered_Min2_samples.tsv", header = T, sep = "\t", fill = T, check.names = F)

  # Save sample IDs
  Samples_Total_before <- colSums(Metaphlan_profiles[,-1])
  table(Samples_Total_before)

  # Rm samples with Total abundance 0
  # The following samples have 0 abundances; Probably only had singletons
  SamplestoRm <- names(which(Samples_Total_before == 0))
  names(which(Samples_Total_before== 0))

  # rm failed samples: Mv this up later
  Metaphlan_profiles <- Metaphlan_profiles %>%
    select(c("taxonomy",!all_of((SamplestoRm)))) # breaks code?
  # if none found then proceed without filtering Samples out

  # subset Metadata df
  Metadata_df_sub <- Metadata_df %>%
    filter(RT_category == "LRT")

# Step 1: Data preparation (Common for all sample types)
  # Get a list of M4 profile samples to retain
  retainSamples <- intersect((colnames(Metaphlan_profiles)[-1]), Metadata_df_sub$SampleID)
    # Check how many are common between both tables
   length(retainSamples)

  # Prepare feature abundances be provided as a feature-by-sample matrix
  row.names(Metaphlan_profiles) <- Metaphlan_profiles$taxonomy
  Metaphlan_profiles <- Metaphlan_profiles[retainSamples]
    # sanity check
    dim(Metaphlan_profiles)

  # subset metadata df to match
  Metadata_df_sub <- Metadata_df_sub %>% filter(SampleID %in% retainSamples)
  row.names(Metadata_df_sub) <- Metadata_df_sub$SampleID
    # sanity check
    dim(Metadata_df_sub)

# Step 2: Numeric values formatting
  # Convert all columns of a data frame to numeric or integer
  Metaphlan_profiles[] <- lapply(Metaphlan_profiles, function(x) {
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
    Sum_abundances <- colSums(Metaphlan_profiles)
    check_feature_abd(Metaphlan_profiles) 
    # will FAIL because MMUPHIN expects a different scale of sample values. The total sum per sample is 100. Therefore, we scale it to 1.

  # Scaling the data from 0-100 range to 0-1
    Metaphlan_profiles_scaled <- Metaphlan_profiles / 100
    check_feature_abd(Metaphlan_profiles_scaled)
    
  # Sanity check batches
    table(Metadata_df_sub$BioProject, Metadata_df_sub$SampleType)
    table(Metadata_df_sub$SampleType)

  # Check if we match sample names between metadata and profiles
  # Step 3.a: Ensure Metadata_df_sub row names match colnames of Metaphlan_profiles_scaled
  Metadata_df_sub <- Metadata_df_sub[match(colnames(Metaphlan_profiles_scaled), row.names(Metadata_df_sub)), ]

  # Step 3.b: Check if order is now identical
  identical(colnames(Metaphlan_profiles_scaled), row.names(Metadata_df_sub))  
  # Should return TRUE

  # add factors
  Metadata_df_sub$BioProject <- factor(Metadata_df_sub$BioProject)

  # sanity check
  table(Metadata_df_sub$BioProject)
  table(Metadata_df_sub$SampleType)
  levels(Metadata_df_sub$BioProject)


# Step 4: Batch correction
  # For single sampleType we don't subset 

  # Run MMUPHin batch correction for this sample type
  Simple_fit_Metaphlan_profiles <- adjust_batch(
    feature_abd = Metaphlan_profiles_scaled,  # Taxa abundance table scaled to 1
    batch = "BioProject",        # Batch effect to correct
    #covariates = c("Healthy"), # Healthy doesn't work because the BioProject and Healthy are colinear i.e, Different BioProject have only Healthy or Diseased samples. 
    data = Metadata_df_sub_, # Metadata file
    control = list(verbose = FALSE)
  )
  
  # Get the corrected abundance table now
  Metaphlan_profiles_scaled_adj <- data.frame(Simple_fit_Metaphlan_profiles$feature_abd_adj)

# Step 5: Save results
  Metaphlan_profiles_scaled_adj$Species <- row.names(Metaphlan_profiles_scaled_adj) 
  write.table(Metaphlan_profiles_scaled_adj, "../MetaPhlan4_results/MiPORT_LRT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", row.names = T, sep = "\t")

# Step 6: Check effect reduction
  # Test variance explained now
  library(vegan, quietly = TRUE)
  set.seed(8284)

  # random check
  Samples_Total_before <- colSums(Metaphlan_profiles_scaled)
  Samples_Total_after <- colSums(Metaphlan_profiles_scaled_adj)

  # calc BC dist
  D_before <- vegdist(t(Metaphlan_profiles_scaled))
  D_after <- vegdist(t(Metaphlan_profiles_scaled_adj))

  # Calculate R-sq with Adonis fit to predict Bio-project
  fit_adonis_before <- adonis2(D_before ~ BioProject, data = Metadata_df_sub, parallel = 8) 
    # For LRT remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

  fit_adonis_after <- adonis2(D_after ~ BioProject, data = Metadata_df_sub, parallel = 8) 
    # For LRT remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

  # Check effect size reduction
  print(fit_adonis_before)
  print(fit_adonis_after)


######################################################################################
##################### Batch correction for URT Sample Type ###########################
######################################################################################


# Read M4 profile file for URT category
Metaphlan_profiles <- read.table("../MetaPhlan4_results/MiPORT_URT_featureTable_species_filtered_Min2_samples.tsv", header = T, sep = "\t", fill = T, check.names = F)

  # NOTE: subset for Buccal_mucosa & Saliva samples since these are the only 2 with multiple batches

  Samples_Total_before <- round(colSums(Metaphlan_profiles[,-1]),1)
  table(Samples_Total_before)

  # Following samples have 0 abundances; Probably only had singletons
  SamplestoRm <- names(which(Samples_Total_before == 0))
  names(which(Samples_Total_before== 0))

  # rm failed samples: Mv this up later
  Metaphlan_profiles <- Metaphlan_profiles %>%
    select(c("taxonomy",!all_of((SamplestoRm)))) # breaks code?
  # if none found then proceed without filtering Samples out

  # Read QC sample ID & metadata list
  Metadata_df <- read.table("../MiPORT_Metadata_downstream_filtered_M4_passed.tsv", header = T, sep = "\t", check.names = F)

  # remove blank rows
  Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]
  Metadata_df <- Metadata_df %>% filter(SampleType != 'Anterior_nares')

  table(Metadata_df$RT_category)

  # subset Metadata df
  Metadata_df_sub <- Metadata_df %>%
    filter(RT_category == "URT")

# Step 1: Data preparation (Common for all sample types)
  # Get a list of M4 profile samples to retain
  retainSamples <- intersect((colnames(Metaphlan_profiles)[-1]), Metadata_df_sub$SampleID)
  # Check how many are common between both tables
  length(retainSamples)

  # Prepare feature abundances be provided as a feature-by-sample matrix
  row.names(Metaphlan_profiles) <- Metaphlan_profiles$taxonomy
  Metaphlan_profiles <- Metaphlan_profiles[retainSamples]

  # sanity check
  dim(Metaphlan_profiles)

  # subset metadata df to match
  Metadata_df_sub <- Metadata_df_sub %>% filter(SampleID %in% retainSamples)
  row.names(Metadata_df_sub) <- Metadata_df_sub$SampleID

  # sanity check
  dim(Metadata_df_sub)

# Step 2: Numeric values formatting
  # Convert all columns of a data frame to numeric or integer
  Metaphlan_profiles[] <- lapply(Metaphlan_profiles, function(x) {
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
    Sum_abundances <- colSums(Metaphlan_profiles)
    check_feature_abd(Metaphlan_profiles) 
    # will FAIL because MMUPHIN expects a different scale of sample values. The total sum per sample is 100. Therefore, we scale it to 1.

  # Scaling the data from 0-100 range to 0-1
    Metaphlan_profiles_scaled <- Metaphlan_profiles / 100
    check_feature_abd(Metaphlan_profiles_scaled)

  # Sanity check batches
    table(Metadata_df_sub$BioProject, Metadata_df_sub$SampleType)
    table(Metadata_df_sub$SampleType)

  # Subset for sample type with multiple batches
  # Keep only these 2 sampleTypes
  SampleTypes_to_Correct <- c("Buccal_mucosa", "Saliva")
  Metadata_df_sub_URT <- Metadata_df_sub %>%
    filter(SampleType %in% c(SampleTypes_to_Correct))

  # sanity check
  table(Metadata_df_sub_URT$SampleType)

  # Get a list of M4 profile samples to retain
  retainSamples <- intersect((colnames(Metaphlan_profiles_scaled)), Metadata_df_sub_URT$SampleID)
  length(retainSamples)
  dim(Metaphlan_profiles_scaled) # before 1497 1630
  
  # subset M4 profiles
  Metaphlan_profiles_scaled_sub <- Metaphlan_profiles_scaled[retainSamples]
  dim(Metaphlan_profiles_scaled_sub) # after 1497 1342 (1027 + 315)

  # Check if we match sample names between metadata and profiles
  # Step 3.a: Ensure Metadata_df_sub_URT row names match colnames of Metaphlan_profiles_scaled_sub
    Metadata_df_sub_URT <- Metadata_df_sub_URT[match(colnames(Metaphlan_profiles_scaled_sub), row.names(Metadata_df_sub_URT)), ]

  # Step 3.b: Check if order is now identical
  identical(colnames(Metaphlan_profiles_scaled_sub), row.names(Metadata_df_sub_URT))  
  # Should return TRUE

  # add factors
  Metadata_df_sub_URT$BioProject <- factor(Metadata_df_sub_URT$BioProject)

  # sanity check
  table(Metadata_df_sub_URT$BioProject)
  table(Metadata_df_sub_URT$SampleType)

  # subset for 
  # Metadata_df_sub <- Metadata_df %>% filter(BioProject != "PRJNA659860")

# Step 4: Batch correction
  # URT has 2 sampletypes to correct. We subset and then correct for each

  # Make a copy of the input URT df
  Metaphlan_profiles_scaled_sub_OG <- Metaphlan_profiles_scaled_sub
  Metadata_df_sub_URT_OG <- Metadata_df_sub_URT

  # A. Correct for Saliva 
    # subset for 1 sampleType in this
    SamplesIDs_toRetain <- Metadata_df_sub_URT_OG %>% 
      filter(SampleType == "Saliva") %>% # Either Buccal_mucosa or Saliva
      pull(SampleID)

  # sanity check
  length(SamplesIDs_toRetain)

  # subset metadata and metaphlan table for these
  Metaphlan_profiles_scaled_sub <- Metaphlan_profiles_scaled_sub_OG[,SamplesIDs_toRetain]
  Metadata_df_sub_URT <- Metadata_df_sub_URT_OG %>% 
    filter(SampleType == "Saliva") 

  # drop unused factors
  Metadata_df_sub_URT$BioProject <- droplevels(Metadata_df_sub_URT$BioProject)

  # sanity check
  dim(Metaphlan_profiles_scaled_sub)
  dim(Metadata_df_sub_URT) 
  table(Metadata_df_sub_URT$BioProject)

  # Run MMUPHin batch correction for this sample type
  Simple_fit_Metaphlan_profiles <- adjust_batch(
    feature_abd = Metaphlan_profiles_scaled_sub,  # Taxa abundance table scaled to 1
    batch = "BioProject",        # Batch effect to correct
    #covariates = c("Healthy"), 
    data = Metadata_df_sub_URT, # Metadata file
    control = list(verbose = FALSE)
  )

  # Get the corrected abundance table for saliva
  Metaphlan_profiles_scaled_adj_Saliva <- data.frame(Simple_fit_Metaphlan_profiles$feature_abd_adj)

  # B. Correct for Saliva 
    # subset for 1 sampleType in this
    SamplesIDs_toRetain <- Metadata_df_sub_URT_OG %>% 
      filter(SampleType == "Buccal_mucosa") %>% # Either Buccal_mucosa or Saliva
      pull(SampleID)

    # sanity check
    length(SamplesIDs_toRetain)

  # subset metadata and metaphlan table for these
  Metaphlan_profiles_scaled_sub <- Metaphlan_profiles_scaled_sub_OG[,SamplesIDs_toRetain]
  Metadata_df_sub_URT <- Metadata_df_sub_URT_OG %>% 
    filter(SampleType == "Buccal_mucosa") 

    # drop unused factors
    Metadata_df_sub_URT$BioProject <- droplevels(Metadata_df_sub_URT$BioProject)

  # sanity check
  dim(Metaphlan_profiles_scaled_sub)
  dim(Metadata_df_sub_URT) 
  table(Metadata_df_sub_URT$BioProject)

  # Run MMUPHin batch correction for this sample type
  Simple_fit_Metaphlan_profiles <- adjust_batch(
    feature_abd = Metaphlan_profiles_scaled_sub,  # Taxa abundance table scaled to 1
    batch = "BioProject",        # Batch effect to correct
    #covariates = c("Healthy"), 
    data = Metadata_df_sub_URT, # Metadata file
    control = list(verbose = FALSE)
  )

  # Get the corrected abundance table now
  Metaphlan_profiles_scaled_adj_Buccal_mucosa <- data.frame(Simple_fit_Metaphlan_profiles$feature_abd_adj)


  # c. Finally, we merge all 3 tables together 
    # Metaphlan_profiles_scaled_adj_Buccal_mucosa,
    # Metaphlan_profiles_scaled_adj_Saliva, 
    # OG_abundance from other sampleTypes (other than Buccal and saliva)

  # quick check
  table(Metadata_df_sub$SampleType)

  # get sampleIds which are not in SampleTypes_to_Correct
  SamplesIDs_toRetain <- Metadata_df_sub %>% 
    filter(!(SampleType %in% SampleTypes_to_Correct)) %>%
    pull(SampleID)

  length(SamplesIDs_toRetain) # Should be 288 (27 + 15 + 246)

  # Merge
  Metaphlan_profiles_scaled_adj <- cbind(Metaphlan_profiles_scaled[, SamplesIDs_toRetain], Metaphlan_profiles_scaled_adj_Buccal_mucosa, Metaphlan_profiles_scaled_adj_Saliva)

  dim(Metaphlan_profiles_scaled_adj)
  dim(Metaphlan_profiles_scaled)

  # Step 5: Save results
  write.table(Metaphlan_profiles_scaled_adj, "../MetaPhlan4_results/MiPORT_URT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", row.names = T, sep = "\t", quote = F)

# Step 6: Check effect reduction

# Test variance explained now
library(vegan, quietly = TRUE)
# Example: Convert "RTMicrobiome" using ASCII values of initials (R=82, T=84)
set.seed(8284)

  # random check
  # A. Buccal mucosa adjustments
    dim(Metaphlan_profiles_scaled_sub)
    dim(Metaphlan_profiles_scaled_adj_Buccal_mucosa)

    Samples_Total_before <- round(colSums(Metaphlan_profiles_scaled), 1)
    Samples_Total_after <- round(colSums(Metaphlan_profiles_scaled_adj_Buccal_mucosa), 1)

    table(Samples_Total_after)
    table(Samples_Total_before)

    # calc BC dist
    D_before <- vegdist(t(Metaphlan_profiles_scaled_sub))
    D_after <- vegdist(t(Metaphlan_profiles_scaled_adj_Buccal_mucosa))

    # Calculate R2 with Adonis fit to predict Bio-project
    fit_adonis_before <- adonis2(D_before ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 
    # Add 'strata = Metadata_df_sub$SampleType' param 

    fit_adonis_after <- adonis2(D_after ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 

    print(fit_adonis_before)
    print(fit_adonis_after)


  # B. Saliva adjustments
    dim(Metaphlan_profiles_scaled_adj_Saliva)
    dim(Metadata_df_sub_URT)
    
    Samples_Total_before <- round(colSums(Metaphlan_profiles_scaled), 1)
    Samples_Total_after <- round(colSums(Metaphlan_profiles_scaled_adj_Saliva), 1)

    table(Samples_Total_after)
    table(Samples_Total_before)

    D_before <- vegdist(t(Metaphlan_profiles_scaled_sub))
    D_after <- vegdist(t(Metaphlan_profiles_scaled_adj_Saliva))

      # Calculate R2 with Adonis fit to predict Bio-project
    fit_adonis_before <- adonis2(D_before ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 
    # Add 'strata = Metadata_df_sub$SampleType' param 

    fit_adonis_after <- adonis2(D_after ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 

    print(fit_adonis_before)
    print(fit_adonis_after)


  # C. Global URT
    dim(Metaphlan_profiles_scaled_sub)
    dim(Metaphlan_profiles_scaled_adj)

    Samples_Total_before <- round(colSums(Metaphlan_profiles_scaled), 1)
    Samples_Total_after <- round(colSums(Metaphlan_profiles_scaled_adj), 1)

    table(Samples_Total_after)
    table(Samples_Total_before)

    # calc BC dist
    D_before <- vegdist(t(Metaphlan_profiles_scaled))
    D_after <- vegdist(t(Metaphlan_profiles_scaled_adj))

    # Calculate R2 with Adonis fit to predict Bio-project
    fit_adonis_before <- adonis2(D_before ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 
    # Add 'strata = Metadata_df_sub$SampleType' param 

    fit_adonis_after <- adonis2(D_after ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 

    print(fit_adonis_before)
    print(fit_adonis_after)

# ---------------------- DONE ---------------------- #


# Side track?: PatientID explains most variance
fit_adonis_before <- adonis2(D_before ~ PatientID, data = Metadata_df_sub_URT, parallel = 8) 
# Remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

fit_adonis_after <- adonis2(D_after ~ PatientID, data = Metadata_df_sub_URT, parallel = 8) # Remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

print(fit_adonis_before)
print(fit_adonis_after)

