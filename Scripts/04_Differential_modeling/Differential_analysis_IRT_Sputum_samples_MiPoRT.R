# Script to detect Microbiol associations with MaAsLin 3
# Takes in batch Corrected LRT samples dataset and metadata table to check for associations. 
# Run LRT model 2 per batch

# Load libraries
for (lib in c('maaslin3', 'dplyr', 'ggplot2', 'knitr', 'tidyr', 'purrr')){
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# The user must supply 
#  a table of per-sample feature abundances (with zeros still included), 
#  a table of per-sample metadata, and 
#  either a formula or fixed effects (random effects may also be added) that specify how the metadata should relate to the feature prevalence (how likely the feature is to be present or absent) and abundance (how much of the feature is there if it's there)

# 1. Read feature table
Metaphlan_profiles <- read.table("../../MetaPhlan4_results/RT_based_feature_table_subsets/MiPORT_IRT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", header = T, sep = "\t")

  # Get number of samples
  nSamples <- (dim(Metaphlan_profiles)[2])
  nTaxa <- dim(Metaphlan_profiles)[1]

  # sanity checks
  print(paste0("Check if these are the number of samples you have: ",nSamples))
  print(paste0("Check if these are the number of taxa you have: ",nTaxa))

  # Read species filter file
  Species_filter_df <- 
    read.table("IRT_Prevalence_min3_sample_species_any_DiseaseGrp.txt",
               sep = '\t',
               header = T)
  
  # Check number of species
  length(unique(Species_filter_df$Species))
    # 276 species are present in min 3 samples per disease group
  
  # Filter metaphlan table to retain species present in at least 10 samples of any group
  
  # filter feature table for these taxa of interest
  Species_filter_names <- unique(Species_filter_df$Species) 
  
  # apply species filter
  Metaphlan_profiles <- Metaphlan_profiles %>%
    filter(row.names(Metaphlan_profiles) %in% Species_filter_names)
    
# 2. Read metadata table and filter for samples if interest
  Metadata_Df <- 
    read.table("../../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", header = T, sep = "\t")

  # remove blank rows from Metadata
  Metadata_Df <- Metadata_Df[rowSums(is.na(Metadata_Df)) < ncol(Metadata_Df), ]

  # Filter metadata df for IRT samples only
  Metadata_Df <- Metadata_Df %>% 
    filter(RT_category == 'IRT')

  # Update rownames with sample IDs
  rownames(Metadata_Df) <- Metadata_Df$Sample

  # Filter Metadata for diseases of interest
  table(Metadata_Df$Disease, Metadata_Df$SampleTypev2, useNA = 'ifany')

'
                      Other Sputum Supraglottal Tongue_dorsum
  Asthma                  0      4            0             0
  COPD                    0      3            0             0
  Covid-19                0      0           59             0
  Cystic Fibrosis         0    162            0             0
  Healthy                49     15            0           418
  Lower lobe collapse     0      1            0             0
  Multiple_diseases       0      4            0             0
  Other                   0      1            0             0
  SMOKING                 0      3            0             0
'

# Make a list of diseases you want to check the association for
  Disease_of_Interest <- c("Healthy", "Cystic Fibrosis")
  Metadata_Df_DOI <- Metadata_Df %>% 
    filter(Disease %in% Disease_of_Interest) %>%
    filter(Richness > 2) %>% # rm samples with observed diversity <2
    filter(SampleTypev2 == "Sputum")

  # Update rownames with sample IDs
  rownames(Metadata_Df_DOI) <- Metadata_Df_DOI$SampleID

  dim(Metadata_Df_DOI)
  #[1] 152  32

  table(Metadata_Df_DOI$Disease, useNA = 'ifany')
  '
  Cystic Fibrosis         Healthy 
              138             14
  '

  # Check sample batches
  table(Metadata_Df_DOI$BioProject, Metadata_Df_DOI$Disease)
    '
              Cystic Fibrosis Healthy
  PRJEB9034                 0      10
  PRJNA316056               5       0
  PRJNA316588               4       4
  PRJNA516442              51       0
  PRJNA516870              72       0
  PRJNA644285               6       0
    ' 

# Create a new grouping variable. Put all Healthy in one group
  Metadata_Df_DOI_sub <- Metadata_Df_DOI %>%
      # Add col for Total read counts (Total reads will be used as a Random effect variable for Maaslin model)
    mutate(Total_reads = After_QC_R1 + After_QC_R2) %>%
    select(all_of(
      c('SampleID', 'BioProject', 'Total_reads', 'PatientID', 'Disease')
      )) %>%
      mutate(Dataset_Group = case_when(
        BioProject == "PRJNA316056" ~ "CF-Dataset-1",
        BioProject == "PRJNA316588" & Disease == "Cystic Fibrosis" ~ "CF-Dataset-2",
        BioProject == "PRJNA516442" ~ "CF-Dataset-3",
        BioProject == "PRJNA516870" ~ "CF-Dataset-4",
        BioProject == "PRJNA644285" ~ "CF-Dataset-5",
        (BioProject == "PRJEB9034" & Disease == "Healthy") ~ "Healthy-Grp",
        (BioProject == "PRJNA316588" & Disease == "Healthy") ~ "Healthy-Grp",
        TRUE ~ "Other"
      ))

    # sanity check
    table(Metadata_Df_DOI_sub$Dataset_Group, useNA = 'ifany')

    # sanity check
    dim(Metadata_Df_DOI_sub)
      # [1] 152  6

# 4. Filter out Metaphlan samples profiles to match metadata DOI sample list
  Metaphlan_profiles_DOI <- Metaphlan_profiles %>% 
    select(all_of(Metadata_Df_DOI_sub$SampleID))

  dim(Metaphlan_profiles_DOI)
  # 152 samples ->  14 Healthy + 138 CF
  # [1] 276 152

    # SKIP [BEGIN]
    # Instead read the output file from this section and filter species
      # 5. Check taxa prevalence and filter taxa, filter samples
      # Count frequency of non-zero species abundances
      Zero_Status_Metaphlan_profiles <- Metaphlan_profiles_DOI == 0
      
      # Convert logical TRUE (absence) and FALSE (presence) to numeric 0 and 1
      Presence_Absence_Df <- as.data.frame(!Zero_Status_Metaphlan_profiles) * 1
      Presence_Absence_Df$Species <- row.names(Presence_Absence_Df)
      
      # Convert wide PA table to long
      PA_long <- Presence_Absence_Df %>%
        pivot_longer(
          cols = -Species,
          names_to = "SampleID",
          values_to = "Presence"
        )
      
      # Join with metadata
      PA_annotated <- PA_long %>%
        left_join(Metadata_Df_DOI, by = "SampleID") %>%
        filter(Disease %in% Disease_of_Interest)
      
      # Calculate prevalence for each species × disease group
      Prevalence_By_Group <- PA_annotated %>%
        group_by(Species, Disease) %>%
        summarise(
          n_present = sum(Presence),
          n_total = n(),
          Prevalence_Percent = round((n_present / n_total) * 100, 2),
          .groups = "drop"
        ) %>%
        arrange(desc(Prevalence_Percent))
      
      # View result
      head(Prevalence_By_Group)
      
      # Filter species present in ≥3 samples in any disease group
      Filtered_Prevalence <- Prevalence_By_Group %>%
        filter(n_present >= 3)
      
      # Optional: Get full list of species passing filter
      Species_to_keep <- Filtered_Prevalence %>%
        distinct(Species) %>%
        pull(Species)
      
      # If you want to apply this filter to your original data:
      Filtered_PA <- Prevalence_By_Group %>%
        filter(Species %in% Species_to_keep)
    
      # save this table
      write.table(Filtered_PA, "IRT_Prevalence_min3_sample_species_any_DiseaseGrp.txt",
                  sep = '\t', row.names = F)
      
        # Preview filtered table
        head(Filtered_PA)
        # Read filter table
        Species_filter <- 
          read.table("IRT_Prevalence_min3_sample_species_any_DiseaseGrp.txt",
                     sep = '\t',
                     header = T)
        
        # Check number of species
        length(unique(Species_filter$Species))
        # 276 species are present in min 3 samples per disease group
        
        # Filter metaphlan table to retain species present in at least 10 samples of any group
        
        # filter feature table for these taxa of interest
        Species_filter_names <-unique(Species_filter$Species) 
        
        dim(Metaphlan_profiles_DOI)
        # filter results: 1497  species -> 276 species found in at least 3 samples of each disease group
        ### other numbers: 103 species found in at least 10 samples of any group
        
      # SKIP [END]
  

# 5: Add factors to metadata groups of interest

  # Loop over each Dataset_Group & run Maaslin3
  Unique_Runs <- unique(Metadata_Df_DOI_sub$Dataset_Group)
  Unique_Runs
  
  levels(Metadata_Df_DOI_sub$Disease)
  unique((Metadata_Df_DOI_sub$Disease))
  
  # add factors to disease
  Metadata_Df_DOI_sub$Disease <-
    factor(Metadata_Df_DOI_sub$Disease,
           levels = c('Healthy', 'Cystic Fibrosis'))
  
  levels(Metadata_Df_DOI_sub$Disease)
  table(Metadata_Df_DOI_sub$Disease)
  
  # Loop
  for (each_Grp in Unique_Runs[-1]) {
    
    # subset metadata
    Metadata_Df_DOI_sub_dataset <- Metadata_Df_DOI_sub %>%
      filter(Dataset_Group %in% c(each_Grp, 'Healthy-Grp'))
    
    # Drop levels and check factors
    Metadata_Df_DOI_sub_dataset$Disease <- 
      droplevels(Metadata_Df_DOI_sub_dataset$Disease)
    
    levels(Metadata_Df_DOI_sub_dataset$Disease)  
    
    # sanity check
    print(table(Metadata_Df_DOI_sub_dataset$Dataset_Group, 
                useNA = 'ifany'))
    
    # Add row names
    row.names(Metadata_Df_DOI_sub_dataset) <- Metadata_Df_DOI_sub_dataset$SampleID
    
    # Filter out Metaphlan samples profiles to match metadata DOI sample list
    Metaphlan_profiles_DOI <- Metaphlan_profiles %>% 
      filter(row.names(Metaphlan_profiles) %in% Species_filter_names) %>%
      select(all_of(Metadata_Df_DOI_sub_dataset$SampleID))

    # Check how the dataframes look like
    Metaphlan_profiles_DOI[1:5, 1:5]
    Metadata_Df_DOI_sub_dataset[1:5, 1:5]
    
    Disease_Name <- levels(Metadata_Df_DOI_sub_dataset$Disease)[2]
    
    # Give a name for directory
    output_Name <- paste0(each_Grp, "-Model-2")
    
    # Fit a model with ~ Disease + BioProject + Total_reads
    RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                           input_metadata = Metadata_Df_DOI_sub_dataset,
                           output = output_Name,
                           formula = '~ Disease + Total_reads',
                           normalization = 'TSS',
                           transform = 'LOG',
                           augment = TRUE,
                           standardize = TRUE,
                           max_significance = 0.1,
                           median_comparison_abundance = TRUE,
                           median_comparison_prevalence = TRUE,
                           max_pngs = 500,
                           cores = 6,
                           coef_plot_vars = c(paste0("Disease ", Disease_Name)),
                           heatmap_vars = c("Total_reads")
                           
    )
    
    }
      
  # SKIP
    # Add factors for BioProject
    Metadata_Df_DOI$BioProject <- factor(Metadata_Df_DOI$BioProject)
      # Update 
    # Add factors for AgeGrp
    table(Metadata_Df_DOI$AgeGroup, useNA = 'ifany')
    '
            Child Young adult       Adult Older adult        <NA> 
             15           7          16         107         147 
      '
    # Replace NA values in AgeGroup with "Unknown"
    Metadata_Df_DOI$AgeGroup[is.na(Metadata_Df_DOI$AgeGroup)] <- "Unknown"
    table(Metadata_Df_DOI$AgeGroup, useNA = 'ifany')
    
    # init vars
    AgeGrp_Factor<- c("Child", "Young_adult", "Adult", "Older_adult", "Unknown")
    AgeGrp_Labels<- c("Child", "Young adult", "Adult", "Older adult", "Unknown")
    
    # add factors
    Metadata_Df_DOI$AgeGroup <- factor(Metadata_Df_DOI$AgeGroup, 
                                       levels = AgeGrp_Factor,
                                       labels = AgeGrp_Labels)
    
    table(Metadata_Df_DOI$AgeGroup, useNA = 'ifany')
    levels(Metadata_Df_DOI$AgeGroup)
  
    
    # set factors for Healthy
    Metadata_Df_DOI$Healthy <- factor(Metadata_Df_DOI$Healthy, 
                                      levels = c("TRUE", "FALSE"),
                                      labels = c("Healthy", "Diseased"))
    
    table(Metadata_Df_DOI$Healthy)
    levels(Metadata_Df_DOI$Healthy)
    
    # factor for Country
    table(Metadata_Df_DOI$Country, Metadata_Df_DOI$Disease, useNA = 'ifany')
    '
           Covid-19 Healthy Pneumonia
  China         0       8        15
  Norway        0      24         0
  UK            0       0       139
  USA         105       0         0
  '
    Metadata_Df_DOI$Country <- factor(Metadata_Df_DOI$Country)
    
    levels(Metadata_Df_DOI$Country)
    
    
    # SKIP

  # factor for disease groups
  table(Metadata_Df_DOI$Disease, useNA = 'ifany')
  
  levels(Metadata_Df_DOI$Disease)
  Metadata_Df_DOI$Disease <- 
    factor(Metadata_Df_DOI$Disease,
           levels = Disease_of_Interest)
  
  table(Metadata_Df_DOI$Disease, useNA = 'ifany')


# Check how the dataframes look like
Metaphlan_profiles_DOI[1:5, 1:5]
Metadata_Df_DOI[1:5, 1:5]


# 7. Fit Maaslin model to associate taxa to Disease

# Model fitting guidelines:
👉 Order does not matter for fixed effects and random effects in a linear mixed model.
👉 All terms are estimated simultaneously, and their contribution is determined by model fitting, not sequence.
👉 If you suspect collinearity or confounding, test different models instead of relying on term order. 🚀

# Model 1 (All) with Patient stratification
# Adding BioProject to reduce batch effect, 
# stratifying with PatientID to reduce repeated measures effect; Doesn't work
# Adding AgeGroup doesn't work

# Read table and filter species from feature table with this info
# save this table
Filtered_PA <- 
  read.table("LRT_Prevalence_min10_sample_species_any_DiseaseGrp.txt",
            sep = '\t', header = T)

# filter feature table for these taxa of interest
Metaphlan_profiles_DOI_filtered <- Metaphlan_profiles_DOI %>% 
  filter( row.names(Metaphlan_profiles_DOI) %in% Species_to_keep) 

dim(Metaphlan_profiles_DOI)

# Total_reads as additional covariate 
set.seed(8284)

# Fit a model with ~ Disease + BioProject + Total_reads
RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                       input_metadata = Metadata_Df_DOI,
                       output = 'T2_LRT_Disease_Batch_RC_output',
                       formula = '~ Disease + BioProject + Total_reads',
                       normalization = 'TSS',
                       transform = 'LOG',
                       augment = TRUE,
                       standardize = TRUE,
                       max_significance = 0.1,
                       median_comparison_abundance = TRUE,
                       median_comparison_prevalence = TRUE,
                       max_pngs = 500,
                       cores = 6,
                       coef_plot_vars = c("Disease Pneumonia",
                                          "Disease Covid-19"),
                       heatmap_vars = c("Total_reads")
                       
)