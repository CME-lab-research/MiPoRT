library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)


# Read M4 profiles; Remove the comment char before reading the file
Metaphlan_profiles <- read.table("../MetaPhlan4_results/Metaphlan4_species_merged_filtered.tsv", header = T, sep = "\t", fill = T, check.names = F)

# Read QC sample ID & metadata list
Metadata_df <- read.table("../MiPORT_Metadata_downstream_filtered_M4_passed.tsv", header = T, sep = "\t", check.names = F)

# repeated measures metadata
#Metadata_df_RM <- read.table("../MiPORT_Metadata_downstream_repeated_measure_filtered.tsv", header = T, sep = "\t")

# remove blank rows
Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]
Metadata_df <- Metadata_df %>% filter(SampleType != 'Anterior_nares')

# Get a list of M4 profile samples to retain
retainSamples <- intersect((colnames(Metaphlan_profiles)[-1]), Metadata_df$SampleID)
#RM_retainSamples <- intersect((colnames(Metaphlan_profiles)[-1]), Metadata_df_RM$SampleID)

# Prepare feature abundances be provided as a feature-by-sample matrix
row.names(Metaphlan_profiles) <- Metaphlan_profiles$taxonomy
Metaphlan_profiles <- Metaphlan_profiles[retainSamples]

# subset metadata also
Metadata_df <- Metadata_df %>% filter(SampleID %in% retainSamples)
row.names(Metadata_df) <- Metadata_df$SampleID

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
  
  # Since the sum is 100 we scale it to 1
  Metaphlan_profiles_fixed <- Metaphlan_profiles / 100
  check_feature_abd(Metaphlan_profiles_fixed)
  
# Check batches
table(Metadata_df$BioProject, Metadata_df$SampleType)
table(Metadata_df$SampleType)

# sanity check to match sample names between metadata and profiles
# Step 1: Ensure Metadata_df row names match colnames of Metaphlan_profiles_fixed
Metadata_df <- Metadata_df[match(colnames(Metaphlan_profiles_fixed), row.names(Metadata_df)), ]

# Step 2: Check if order is now identical
identical(colnames(Metaphlan_profiles_fixed), row.names(Metadata_df))  # Should return TRUE

# add factors
Metadata_df$BioProject <- factor(Metadata_df$BioProject)

# subset
Metadata_df_sub <- Metadata_df %>% filter(BioProject !="PRJNA659860")
Metadata_df_sub$BioProject <- droplevels(Metadata_df_sub$BioProject)
levels(Metadata_df_sub$BioProject)

Metaphlan_profiles_fixed_sub <- Metaphlan_profiles_fixed[,Metadata_df_sub$SampleID] 

# Adjust for batch effects while controlling for the effect of SampleType
Simple_fit_Metaphlan_profiles <- adjust_batch(
  feature_abd = Metaphlan_profiles_fixed_sub,  # Taxa abundance table scaled to 1
  batch = "BioProject",        # Batch effect to correct
  data = Metadata_df_sub, # Metadata file
  control = list(verbose = FALSE)
)

# Get the corrected abundance table now
Metaphlan_profiles_fixed_adj <- data.frame(Simple_fit_Metaphlan_profiles$feature_abd_adj)

Metaphlan_profiles_fixed_adj$Species <- row.names(Metaphlan_profiles_fixed_adj) 
write.table(Metaphlan_profiles_fixed_adj, "../MetaPhlan4_results/Metaphlan4_species_merged_filtered_batchCorrected.tsv", row.names = T, sep = "\t")

# Test variance explained now
library(vegan, quietly = TRUE)
# Example: Convert "RTMicrobiome" using ASCII values of initials (R=82, T=84)
set.seed(8284)

# calc BC dist
D_before <- vegdist(t(Metaphlan_profiles_fixed_sub))
D_after <- vegdist(t(Metaphlan_profiles_fixed_adj))

# Adonis on Bioproject effect
fit_adonis_before <- adonis2(D_before ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType)

fit_adonis_after <- adonis2(D_after ~ BioProject, data = Metadata_df_sub, parallel = 8, strata = Metadata_df_sub$SampleType)

print(fit_adonis_before)
print(fit_adonis_after)

# Adonis on SampleType effect
fit_adonis_before <- adonis2(D_before ~ SampleType, data = Metadata_df_sub, parallel = 8)

fit_adonis_after <- adonis2(D_after ~ SampleType, data = Metadata_df_sub, parallel = 8)

print(fit_adonis_before)
print(fit_adonis_after)

