library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)


# Read M4 profile file for LRT category
Metaphlan_profiles <- read.table("../MetaPhlan4_results/MiPORT_IRT_featureTable_species_filtered_Min2_samples.tsv", header = T, sep = "\t", fill = T, check.names = F)


Samples_Total_before <- colSums(Metaphlan_profiles[,-1])
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

# repeated measures metadata
#Metadata_df_RM <- read.table("../MiPORT_Metadata_downstream_repeated_measure_filtered.tsv", header = T, sep = "\t")

# remove blank rows
Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]
Metadata_df <- Metadata_df %>% filter(SampleType != 'Anterior_nares')

# subset Metadata df
Metadata_df_sub <- Metadata_df %>%
  filter(RT_category == "IRT")
  
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

# MMUPHin fn
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

# Check if cols are numeric values
  Sum_abundances <- colSums(Metaphlan_profiles)
  check_feature_abd(Metaphlan_profiles) 
  # will FAIL ..,
  # since the sum is 100. Therefore, we scale it to 1
  Metaphlan_profiles_fixed <- Metaphlan_profiles / 100
  check_feature_abd(Metaphlan_profiles_fixed)
  
# Check batches
table(Metadata_df_sub$BioProject, Metadata_df_sub$SampleType)
table(Metadata_df_sub$SampleType)

# sanity check to match sample names between metadata and profiles
# Step 1: Ensure Metadata_df_sub row names match colnames of Metaphlan_profiles_fixed
Metadata_df_sub <- Metadata_df_sub[match(colnames(Metaphlan_profiles_fixed), row.names(Metadata_df_sub)), ]

# Step 2: Check if order is now identical
identical(colnames(Metaphlan_profiles_fixed), row.names(Metadata_df_sub))  # Should return TRUE

# add factors
Metadata_df_sub$BioProject <- factor(Metadata_df_sub$BioProject)

table(Metadata_df_sub$BioProject)
table(Metadata_df_sub$SampleType)

# subset
#Metadata_df_sub <- Metadata_df %>% filter(BioProject !="PRJNA659860")
#Metadata_df_sub$BioProject <- droplevels(Metadata_df_sub$BioProject)
levels(Metadata_df_sub$BioProject)

Metaphlan_profiles_fixed_sub <- Metaphlan_profiles_fixed[,Metadata_df_sub$SampleID] 

# Adjust for batch effects while controlling for the effect of SampleType
Simple_fit_Metaphlan_profiles <- adjust_batch(
  feature_abd = Metaphlan_profiles_fixed,  # Taxa abundance table scaled to 1
  batch = "BioProject",        # Batch effect to correct
  covariates = c("Healthy"),
  data = Metadata_df_sub, # Metadata file
  control = list(verbose = FALSE)
)

$$ For IRT
# We subset for Sputum samples and correct only for that since no other samples types are present in multiple BioProjects
# Select sampleType
sample = "Sputum"

# Subset metadata and feature table for the current sample type
  subset_metadata <- Metadata_df_sub[Metadata_df_sub$SampleType == sample, ]
  subset_features <- Metaphlan_profiles_fixed[, subset_metadata$SampleID]
  
  # Drop levels
  subset_metadata$BioProject <- droplevels(subset_metadata$BioProject)
  
  # Run MMUPHin batch correction for this sample type
  corrected_data <- adjust_batch(
    feature_abd = subset_features,
    batch = "BioProject",           # Correct for dataset effect
    #covariates = c("Healthy"),       # Include health status
    data = subset_metadata,
    control = list(verbose = FALSE)
  )
  
  # Get the corrected abundance table now
  Metaphlan_profiles_fixed_adj <- data.frame(corrected_data$feature_abd_adj)
  
  # Merge new Sputum sample abundances with other sampleType abundances
  # get sampleIds which are not Sputum
  SamplesIDs_toRetain <- Metadata_df_sub %>% 
    filter(SampleType != "Sputum") %>%
    pull(SampleID)
  
  Metaphlan_profiles_fixed_adj <- cbind(Metaphlan_profiles_fixed[, SamplesIDs_toRetain], data.frame(corrected_data$feature_abd_adj))
  
  # Change rownames to colnames
  #Metaphlan_profiles_fixed_adj$Species <- row.names(Metaphlan_profiles_fixed_adj) 
  write.table(Metaphlan_profiles_fixed_adj, "../MetaPhlan4_results/MiPORT_IRT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", row.names = T, sep = "\t", quote = F)
  
$$

# For single sampleType > Get the corrected abundance table now
Metaphlan_profiles_fixed_adj <- data.frame(Simple_fit_Metaphlan_profiles$feature_abd_adj)

Metaphlan_profiles_fixed_adj$Species <- row.names(Metaphlan_profiles_fixed_adj) 
write.table(Metaphlan_profiles_fixed_adj, "../MetaPhlan4_results/MiPORT_LRT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", row.names = T, sep = "\t")

# Test variance explained now
library(vegan, quietly = TRUE)
# Example: Convert "RTMicrobiome" using ASCII values of initials (R=82, T=84)
set.seed(8284)

# random check
Samples_Total_before <- colSums(Metaphlan_profiles_fixed)
Samples_Total_after <- colSums(Metaphlan_profiles_fixed_adj)

# calc BC dist
D_before <- vegdist(t(Metaphlan_profiles_fixed))
D_after <- vegdist(t(Metaphlan_profiles_fixed_adj))

# Calculate R2 with Adonis fit to predict Bio-project
fit_adonis_before <- adonis2(D_before ~ BioProject, data = Metadata_df_sub, parallel = 8) 
  # For LRT remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

fit_adonis_after <- adonis2(D_after ~ BioProject, data = Metadata_df_sub, parallel = 8) 
  # For LRT remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

print(fit_adonis_before)
print(fit_adonis_after)


###### DONE: Repeat for the other RT categories ######

#### For URT

# Read M4 profile file for LRT category
Metaphlan_profiles <- read.table("../MetaPhlan4_results/MiPORT_URT_featureTable_species_filtered_Min2_samples.tsv", header = T, sep = "\t", fill = T, check.names = F)


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

# repeated measures metadata
#Metadata_df_RM <- read.table("../MiPORT_Metadata_downstream_repeated_measure_filtered.tsv", header = T, sep = "\t")

# remove blank rows
Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]
Metadata_df <- Metadata_df %>% filter(SampleType != 'Anterior_nares')

table(Metadata_df$RT_category)

# subset Metadata df
Metadata_df_sub <- Metadata_df %>%
  filter(RT_category == "URT")

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

# MMUPHin fn
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

# Check if cols are numeric values
Sum_abundances <- colSums(Metaphlan_profiles)
check_feature_abd(Metaphlan_profiles) 
# will FAIL ..,
# since the sum is 100. Therefore, we scale it to 1
Metaphlan_profiles_fixed <- Metaphlan_profiles / 100
check_feature_abd(Metaphlan_profiles_fixed)

# Check batches
table(Metadata_df_sub$BioProject, Metadata_df_sub$SampleType)
table(Metadata_df_sub$SampleType)

# subset for Buccal_mucosa & Saliva samples since these are the only 2 with multiple batches

SampleTypes_to_Correct <- c("Buccal_mucosa", "Saliva")
Metadata_df_sub_tmp <- Metadata_df_sub %>%
  filter(SampleType %in% c(SampleTypes_to_Correct))

table(Metadata_df_sub_tmp$SampleType)

# Get a list of M4 profile samples to retain
retainSamples <- intersect((colnames(Metaphlan_profiles_fixed)), Metadata_df_sub_tmp$SampleID)

length(retainSamples)

# Keep only these 2 sampleTypes
dim(Metaphlan_profiles_fixed) # before 1497 1630
Metaphlan_profiles_fixed_sub <- Metaphlan_profiles_fixed[retainSamples]

dim(Metaphlan_profiles_fixed_sub) # before 1497 1342 (1027 + 315)

# sanity check to match sample names between metadata and profiles
# Step 1: Ensure Metadata_df_sub row names match colnames of Metaphlan_profiles_fixed
Metadata_df_sub_tmp <- Metadata_df_sub_tmp[match(colnames(Metaphlan_profiles_fixed_sub), row.names(Metadata_df_sub_tmp)), ]

# Step 2: Check if order is now identical
identical(colnames(Metaphlan_profiles_fixed_sub), row.names(Metadata_df_sub_tmp))  # Should return TRUE

# add factors
Metadata_df_sub_tmp$BioProject <- factor(Metadata_df_sub_tmp$BioProject)

table(Metadata_df_sub_tmp$BioProject)
table(Metadata_df_sub_tmp$SampleType)

# subset for 
#Metadata_df_sub <- Metadata_df %>% filter(BioProject !="PRJNA659860")
#Metadata_df_sub$BioProject <- droplevels(Metadata_df_sub$BioProject)
levels(Metadata_df_sub_tmp$BioProject)

# Has both sampletypes to correct
Metaphlan_profiles_fixed_sub_OG <- Metaphlan_profiles_fixed_sub
Metadata_df_sub_tmp_OG <- Metadata_df_sub_tmp

#SampleTypes_to_Correct <- c("Buccal_mucosa", "Saliva")
#Metadata_df_sub_tmp <- Metadata_df_sub %>%
#  filter(SampleType %in% c(SampleTypes_to_Correct))

# subset for 1 sampleType in this
SamplesIDs_toRetain <- Metadata_df_sub_tmp_OG %>% 
  filter(SampleType == "Saliva") %>% # Either Buccal_mucosa or Saliva
  pull(SampleID)

# sanity check
length(SamplesIDs_toRetain)

# subset metadata and metaphlan table for these
Metaphlan_profiles_fixed_sub <- Metaphlan_profiles_fixed_sub_OG[,SamplesIDs_toRetain]
Metadata_df_sub_tmp <- Metadata_df_sub_tmp_OG %>% 
  filter(SampleType == "Saliva") 

Metadata_df_sub_tmp$BioProject <- droplevels(Metadata_df_sub_tmp$BioProject)

dim(Metaphlan_profiles_fixed_sub)
dim(Metadata_df_sub_tmp) 

table(Metadata_df_sub_tmp$BioProject)

# Adjust for batch effects while controlling for the effect of SampleType
Simple_fit_Metaphlan_profiles <- adjust_batch(
  feature_abd = Metaphlan_profiles_fixed_sub,  # Taxa abundance table scaled to 1
  batch = "BioProject",        # Batch effect to correct
  #covariates = c("Healthy"),
  data = Metadata_df_sub_tmp, # Metadata file
  control = list(verbose = FALSE)
)

# Get the corrected abundance table now
Metaphlan_profiles_fixed_adj_Buccal_mucosa <- data.frame(Simple_fit_Metaphlan_profiles$feature_abd_adj)

# repeat the same for saliva and get
Metaphlan_profiles_fixed_adj_Saliva <- data.frame(Simple_fit_Metaphlan_profiles$feature_abd_adj)

# Finally merge all 3 tables together (Metaphlan_profiles_fixed_adj_Buccal_mucosa,
# Metaphlan_profiles_fixed_adj_Saliva, 
# OG_abundance from other sampleTypes)

table(Metadata_df_sub$SampleType)

# get sampleIds which are not in SampleTypes_to_Correct
SamplesIDs_toRetain <- Metadata_df_sub %>% 
  filter(!(SampleType %in% SampleTypes_to_Correct)) %>%
  pull(SampleID)

length(SamplesIDs_toRetain) # Should be 288 (27 + 15 + 246)

# Merge
Metaphlan_profiles_fixed_adj <- cbind(Metaphlan_profiles_fixed[, SamplesIDs_toRetain], Metaphlan_profiles_fixed_adj_Buccal_mucosa, Metaphlan_profiles_fixed_adj_Saliva)

dim(Metaphlan_profiles_fixed_adj)
dim(Metaphlan_profiles_fixed)

# Change rownames to colnames
#Metaphlan_profiles_fixed_adj$Species <- row.names(Metaphlan_profiles_fixed_adj) 
write.table(Metaphlan_profiles_fixed_adj, "../MetaPhlan4_results/MiPORT_URT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", row.names = T, sep = "\t", quote = F)

# Test variance explained now
library(vegan, quietly = TRUE)
# Example: Convert "RTMicrobiome" using ASCII values of initials (R=82, T=84)
set.seed(8284)

# random check
dim(Metaphlan_profiles_fixed_sub)
#dim(Metaphlan_profiles_fixed_adj_Buccal_mucosa)
dim(Metaphlan_profiles_fixed_adj_Saliva)
dim(Metadata_df_sub_tmp)

dim(Metaphlan_profiles_fixed_adj)
dim(Metadata_df_sub)

Samples_Total_before <- round(colSums(Metaphlan_profiles_fixed), 1)
Samples_Total_after <- round(colSums(Metaphlan_profiles_fixed_adj), 1)

table(Samples_Total_after)
table(Samples_Total_before)

# calc BC dist
D_before <- vegdist(t(Metaphlan_profiles_fixed))
D_after <- vegdist(t(Metaphlan_profiles_fixed_adj))

#D_after <- vegdist(t(Metaphlan_profiles_fixed_adj_Buccal_mucosa))
#D_after <- vegdist(t(Metaphlan_profiles_fixed_adj_Saliva))

# Calculate R2 with Adonis fit to predict Bio-project
fit_adonis_before <- adonis2(D_before ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 
# Add 'strata = Metadata_df_sub$SampleType' param 

fit_adonis_after <- adonis2(D_after ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType) 


print(fit_adonis_before)
print(fit_adonis_after)

# Side track?: PatientID explains most variance
fit_adonis_before <- adonis2(D_before ~ PatientID, data = Metadata_df_sub_tmp, parallel = 8) 
# Remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

fit_adonis_after <- adonis2(D_after ~ PatientID, data = Metadata_df_sub_tmp, parallel = 8) # Remove 'strata = Metadata_df_sub$SampleType' param because LRT has only BAL

print(fit_adonis_before)
print(fit_adonis_after)

