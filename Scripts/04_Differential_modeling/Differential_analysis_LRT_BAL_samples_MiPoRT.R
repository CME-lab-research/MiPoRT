# Script to detect Microbiome association with MaAsLin 3
# Takes in batchCorrected LRT samples dataset and metadata table to check for associations

# Load libraries
for (lib in c('maaslin3', 'dplyr', 'ggplot2', 'knitr', 'tidyr', 'purrr')){
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
  }

# The user must supply 
#  a table of per-sample feature abundances (with zeros still included), 
#  a table of per-sample metadata, and 
#  either a formula or fixed effects (random effects may also be added) that specify how the metadata should relate to the feature prevalence (how likely the feature is to be present or absent) and abundance (how much of the feature is there if it's there)

# 1. Read feature table
Metaphlan_profiles <- read.table("../../MetaPhlan4_results/RT_based_feature_table_subsets/MiPORT_LRT_featureTable_sp_filt_Min2_samples_batchCorrected.tsv", header = T, sep = "\t")
  
# Get number of samples
nSamples <- (dim(Metaphlan_profiles)[2])
nTaxa <- dim(Metaphlan_profiles)[1]
  # sanity checks
  print(paste0("Check if these are the number of samples you have: ",nSamples))
  print(paste0("Check if these are the number of taxa you have: ",nTaxa))

# 2. Read metadata table and filter for samples if interest
Metadata_Df <- read.table("../../MiPORT_Metadata_downstream_filtered_M4_passed_v2.tsv", header = T, sep = "\t")

  # remove blank rows from Metadata
  Metadata_Df <- Metadata_Df[rowSums(is.na(Metadata_Df)) < ncol(Metadata_Df), ]
  
  # Filter metadata df for LRT samples only
  Metadata_Df <- Metadata_Df %>% filter(RT_category == 'LRT')
    # 585 * 32 samples
  
  # Update rownames with sample IDs
  rownames(Metadata_Df) <- Metadata_Df$Sample
  
  # Filter Metadata for diseases of interest
  table(Metadata_Df$Disease, Metadata_Df$SampleTypev2, useNA = 'ifany')
  
  '
                  BAL
  Asthma            4
  COPD             45
  Covid-19        178
  DIP               1
  Healthy          52
  HP                1
  IPF              10
  LTX_BO            8
  LTX_non-BO       12
  Lung cancer       6
  Lung_cancer      23
  Other             1
  Pneumonia       194
  Sarcoidosis      23
  Tuberculosis     22
  Unspecified ILS   5
  '
  
  # Make a list of diseases you want to check the association for
  Disease_of_Interest <- c("Healthy", "Pneumonia", "Covid-19")
  Metadata_Df_DOI <- Metadata_Df %>% 
    filter(Disease %in% Disease_of_Interest) %>%
    filter(observed > 2) # rm samples with observed diversity <2
  
  
  # Update rownames with sample IDs
  rownames(Metadata_Df_DOI) <- Metadata_Df_DOI$SampleID
  
  dim(Metadata_Df_DOI)
  #[1] 292  33
  
  table(Metadata_Df_DOI$Disease, useNA = 'ifany')
    ' N = 
    Covid-19   Healthy Pneumonia 
          105        33       154 
    '

# Check sample batches
table(Metadata_Df_DOI$BioProject, Metadata_Df_DOI$Disease)

'
              Covid-19 Healthy Pneumonia
  PRJEB29011         0       0       139
  PRJNA413615        0       0        15
  PRJNA655567        0       1         0
  PRJNA687506      105       0         0
  PRJNA757846        0       8         0
  UniBergen          0      24         0
' 
# Remove project PRJNA655567 as there is only 1 sample
Metadata_Df_DOI <- Metadata_Df_DOI %>%
  filter(BioProject != 'PRJNA655567')

# Add factors for BioProject  
  # update the Bio-Project factors to give a base line for comparing abundances
  Metadata_Df_DOI %>%
    group_by(BioProject, Healthy) %>%
    summarise(
      Summary_obs = mean(observed),
      Summary_shannon = mean(diversity_shannon),
      Total_Reads = mean(After_QC_R1 + After_QC_R2)
    )
  '
  `summarise()` has grouped output by 'BioProject'. You can override using the
  `.groups` argument.
  # A tibble: 5 × 5
      # Groups:   BioProject [5]
        BioProject  Healthy Summary_obs Summary_shannon Total_Reads
        <chr>       <chr>         <dbl>           <dbl>       <dbl>
      1 PRJEB29011  FALSE          14.3           1.30     4955994.
      2 PRJNA413615 FALSE          15.1           0.609    6329962.
      3 PRJNA687506 FALSE          12.6           1.36    59555373.
      4 PRJNA757846 TRUE           74.5           2.81     5684356.
      5 UniBergen   TRUE           22.6           2.10     1772529.
        
        Since, BioProject PRJNA757846 has the higher observed species counts and shannon diversity (on average) with more read counts we select that project as Healthy baseline.
  '
  levels(Metadata_Df_DOI$BioProject)
  
  Metadata_Df_DOI$BioProject <- factor(Metadata_Df_DOI$BioProject,
                                       levels = c("PRJNA757846", 
                                                  "UniBergen",
                                                  "PRJEB29011",
                                                  "PRJNA413615",
                                                  "PRJNA687506"
                                                  ))
  levels(Metadata_Df_DOI$BioProject)  
    '[1] "PRJNA757846" "UniBergen"   "PRJEB29011"  "PRJNA413615" "PRJNA687506" '
  
  # sanity check
  dim(Metadata_Df_DOI)
  # [1] 291  33

# 3. Add a new metadata variable  
  # Add col for Total read counts (Total reads will be used as a Random effect variable for Maaslin model)
  Metadata_Df_DOI$Total_reads <- Metadata_Df_DOI$After_QC_R1 + Metadata_Df_DOI$After_QC_R2

# 4. Filter out Metaphlan samples profiles to match metadata DOI sample list
  Metaphlan_profiles_DOI <- Metaphlan_profiles %>% 
    select(all_of(Metadata_Df_DOI$SampleID))
  
  dim(Metaphlan_profiles_DOI)
  # 291 samples ->  105 + 32 + 154
    # [1] 1497  291
  
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
  
  # Get full list of species passing filter
  Species_to_keep <- Filtered_Prevalence %>%
    distinct(Species) %>%
    pull(Species)
  
  # If you want to apply this filter to your original data:
  Filtered_PA <- Prevalence_By_Group %>%
    filter(Species %in% Species_to_keep)
  
  # save this table
  write.table(Filtered_PA, "LRT_Prevalence_min3_sample_species_any_DiseaseGrp.txt",
              sep = '\t', row.names = F)
  
  # Preview filtered table
  head(Filtered_PA)
  
  # Filter metaphlan table to retain species present in at least 10 samples of any group
  
  # filter feature table for these taxa of interest
  Metaphlan_profiles_DOI_filtered <- Metaphlan_profiles_DOI %>% 
    filter( row.names(Metaphlan_profiles_DOI) %in% Species_to_keep) 
  
  dim(Metaphlan_profiles_DOI)
  # filter results: 1497  species 
  dim(Metaphlan_profiles_DOI_filtered)
  # filter results: 1497  species 
  #  -> 312 species found in at least 3 samples of any group
  #  -> 103 species found in at least 10 samples of any group
  
  # save this table for quick Maaslin analysis
  write.table(Metaphlan_profiles_DOI_filtered, 
              "LRT_Filtered_sp_Prevalence_min2_perGrp.txt",
              row.names = T, 
              sep = '\t', quote = F)
  
  
# 6. Add factors to Metadata groups
  # Add row names
  row.names(Metadata_Df_DOI) <- Metadata_Df_DOI$SampleID
  
  # check factors for BioProject
  Metadata_Df_DOI$BioProject <- factor(Metadata_Df_DOI$BioProject)
  levels(Metadata_Df_DOI$BioProject)
  
  # SKIP: Add factors for AgeGrp
  table(Metadata_Df_DOI$AgeGroup, useNA = 'ifany')
  '
        Child Young adult       Adult Older adult        <NA> 
         15           7          16         107         147 
  '
    # Replace NA values in AgeGroup with "Unknown"
    Metadata_Df_DOI$AgeGroup[is.na(Metadata_Df_DOI$AgeGroup)] <- "Unknown"
    table(Metadata_Df_DOI$AgeGroup, useNA = 'ifany')
  
    # Skip: init vars
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
  
  # factor for disease groups
  table(Metadata_Df_DOI$Disease, useNA = 'ifany')
  
  levels(Metadata_Df_DOI$Disease)
  Disease_of_Interest
  
  Metadata_Df_DOI$Disease <- factor(Metadata_Df_DOI$Disease,
                                    levels = Disease_of_Interest
                                    )
  
  table(Metadata_Df_DOI$Disease, useNA = 'ifany')
  
  # calculate total Reads
  #Metadata_Df_DOI$Total_reads <- (Metadata_Df_DOI$After_QC_R1 + Metadata_Df_DOI$After_QC_R2)
    
  # Check how the dataframes look like
  Metaphlan_profiles_DOI[1:5, 1:5]
  Metadata_Df_DOI[1:5, 1:5]
  

# 7. Fit Maaslin model to associate taxa to Disease
  
  # Model fitting guidelines:
  👉 Order does not matter for fixed effects and random effects in a linear mixed model.
  👉 All terms are estimated simultaneously, and their contribution is determined by model fitting, not sequence.
  👉 If you suspect collinearity or confounding, test different models instead of relying on term order. 🚀
  
  # Model 1 (All) with Patient stratification
  # Adding BioProject doesn't make sense since it's co linear to Disease groups. We also reduce sample size
  # stratifying with PatientID to reduce repeated measures effect; Doesn't work
  # Adding AgeGroup doesn't work
  
  # Total_reads as additional covariate 
  set.seed(8284)
  
  # Fit a model with ~ Disease + BioProject + Total_reads
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI_filtered,
                      input_metadata = Metadata_Df_DOI,
                      output = 'Final_LRT_Disease_RC_output',
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
                      coef_plot_vars = c("Disease Pneumonia",
                                         "Disease Covid-19"),
                      heatmap_vars = c("Total_reads")
                      )

  # AgeGrp gives errors for all associations so we remove it 

  # With Age group we get errors for all associations: 
  "INFO::All associations had errors or were insignificant."
  
  
  # Model 2 (M2) without Patient stratifications
  RT_fit_out_2 <- maaslin3(input_data = Metaphlan_profiles_DOI_filtered,
                         input_metadata = Metadata_Df_DOI,
                         output = 'Final_M2_Batch_perDisease_output',
                         formula = '~ BioProject + (1 | Disease) ',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = TRUE,
                         max_pngs = 250,
                         cores = 6,
                         heatmap_vars = c("Total_reads")
                         )

  # Model 3 (M3) Disease only
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                           input_metadata = Metadata_Df_DOI,
                           output = 'LRT_M3_output',
                           formula = '~ Disease',
                           normalization = 'TSS',
                           transform = 'LOG',
                           augment = TRUE,
                           standardize = TRUE,
                           max_significance = 0.1,
                           median_comparison_abundance = TRUE,
                           median_comparison_prevalence = TRUE,
                           max_pngs = 250,
                           cores = 4)
  
  # Model 4 (M4) Disease + Healthy + Country + Total_reads
  RT_fit_out_2 <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M4_output',
                         formula = '~ Disease + Healthy + Country + Total_reads',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = TRUE,
                         max_pngs = 250,
                         cores = 4)
  # Bioproject and Disease are confounded because
  # Check sample batches
  table(Metadata_Df_DOI$BioProject, Metadata_Df_DOI$Disease)
  '
                Healthy Pneumonia
  PRJEB29011        0       178
  PRJNA413615       0        16
  PRJNA655567       5         0
  PRJNA757846      14         0
  UniBergen        33         0
  
    '
  # Therefore, any taxa associated to PRJNA413615 & PRJEB29011 could also be due to disease. 
  
  # Model 5 (M5) Disease + AgeGroup + PatientID
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                           input_metadata = Metadata_Df_DOI,
                           output = 'LRT_M5_Dis_Pls_AgeGrp_Pls_PID_output',
                           formula = '~ Disease + AgeGroup + (1 | PatientID)',
                           normalization = 'TSS',
                           transform = 'LOG',
                           augment = TRUE,
                           standardize = TRUE,
                           max_significance = 0.1,
                           median_comparison_abundance = TRUE,
                           median_comparison_prevalence = FALSE,
                           max_pngs = 250,
                           cores = 8)
  
  # Give errors for all associations: 
  "INFO::All associations had errors or were insignificant."
  
  # Model 6 (M6) Disease + AgeGroup 
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M6_Dis_Pls_AgeGrp_output',
                         formula = '~ Disease + AgeGroup',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = FALSE,
                         max_pngs = 250,
                         cores = 8) 
  "INFO::All associations had errors or were insignificant."
  
  table(Metadata_Df_DOI$AgeGroup, Metadata_Df_DOI$Disease)

  # Model 7 (M7) Disease + Gender + Country
  table(Metadata_Df_DOI$Disease, Metadata_Df_DOI$Country, useNA = 'ifany')
  table(Metadata_Df_DOI$Disease, Metadata_Df_DOI$Gender, useNA = 'ifany')
  
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M7_Dis_Pls_Gender_Pls_Country_output',
                         formula = '~ Disease + Gender + Country',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = FALSE,
                         max_pngs = 250,
                         cores = 8) 
  
  "INFO::All associations had errors or were insignificant."
  
  # Model 8 (M8) Disease + Country
  table(Metadata_Df_DOI$Disease, Metadata_Df_DOI$Country, useNA = 'ifany')
  
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M8_Dis_Pls_Country_output',
                         formula = '~ Disease + Country',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = FALSE,
                         max_pngs = 250,
                         cores = 8) 
  
  
  # Model 9 (M9) Country
  table(Metadata_Df_DOI$BioProject, Metadata_Df_DOI$Country, useNA = 'ifany')
  
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M9_Dis_Pls_BioProj_Pls_Country_output',
                         formula = '~ BioProject + Country',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = FALSE,
                         max_pngs = 250,
                         cores = 8) 
  "INFO::All associations had errors or were insignificant."
  
  # Model 10 (M10) Effect of disease with Country stratification
  table(Metadata_Df_DOI$Disease, Metadata_Df_DOI$Country, useNA = 'ifany')
  
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M10_Dis_Pls_BioProj_Pls_Country_output',
                         formula = '~ Disease + BioProject + Country',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = FALSE,
                         max_pngs = 250,
                         cores = 8) 
  
  # Model 11 (M11) Effect of disease with read count
  table(Metadata_Df_DOI$Disease, Metadata_Df_DOI$Country, useNA = 'ifany')
  
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M11_Dis_Pls_Total_reads_output',
                         formula = '~ Disease + Total_reads',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = FALSE,
                         max_pngs = 250,
                         cores = 8) 
  
  # final cov to add: Country, BioProject, Total_reads
  # Model 12 (M12) Effect of disease with read count
  table(Metadata_Df_DOI$Disease, Metadata_Df_DOI$Country, useNA = 'ifany')
  
  RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                         input_metadata = Metadata_Df_DOI,
                         output = 'LRT_M12_Dis_Pls_All_ImpCov_output',
                         formula = '~ Disease + BioProject + Country + Total_reads',
                         normalization = 'TSS',
                         transform = 'LOG',
                         augment = TRUE,
                         standardize = TRUE,
                         max_significance = 0.1,
                         median_comparison_abundance = TRUE,
                         median_comparison_prevalence = FALSE,
                         max_pngs = 250,
                         cores = 8) 

  # Fix this model: Change the BioProject names to Healthy_Dt_1, Healthy_Dt_2, Pneumonia_Dt_1 etc.

  # Define the dataset mapping
    dataset_labels <- c("PRJEB29011" = "Pneumonia_D1",
                        "PRJNA413615" = "Pneumonia_D2",
                        "PRJNA655567" = "Healthy_D1",
                        "PRJNA757846" = "Healthy_D2",
                        "UniBergen" = "Healthy_D3")
  
  # Create a new column in Metadata_df
    Metadata_Df_DOI$DatasetID <- dataset_labels[Metadata_Df_DOI$BioProject]
  
    # add factors
    Metadata_Df_DOI$DatasetID <- factor(Metadata_Df_DOI$DatasetID)
    
    # Use this col as replacement for BioProject
    RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_DOI,
                           input_metadata = Metadata_Df_DOI,
                           output = 'tmpLRT_M12_Dis_Pls_All_ImpCov_output',
                           formula = '~ (1 + Disease) + DatasetID + Country + Total_reads',
                           normalization = 'TSS',
                           transform = 'LOG',
                           augment = TRUE,
                           standardize = TRUE,
                           max_significance = 0.1,
                           median_comparison_abundance = TRUE,
                           median_comparison_prevalence = FALSE,
                           max_pngs = 250,
                           cores = 8) 
    
    # Global BAL sample model: Covs to add > Country, BioProject, Total_reads
    # Model 13 (M13) Effect on Health + Read count
    # Use DatasetID col as replacement for BioProject
    # subset Metaphlan_profiles with filtered taxa # Process getting killed, try filtering more taxa.
    # remove samples with observed diversity > 2
    Samples_filter_2 <- Metadata_Df %>% filter(observed > 2) %>% pull(SampleID)
    
    Metaphlan_profiles_sub <- Metaphlan_profiles %>%
      filter( row.names(Metaphlan_profiles) %in% Taxa_filter) %>% # Filter Taxa present in >2 samples
      select(Samples_filter_2) %>% # Filter samples with obs > 2
      t() %>% data.frame()
    
    # sanity check
    dim(Metaphlan_profiles_sub)
    
    # subset metadata also
    Metadata_Df_sub <- Metadata_Df %>% filter(observed > 2)
    dim(Metadata_Df_sub)
    
    Metadata_Df_sub$Country <- factor(Metadata_Df_sub$Country)
    Metadata_Df_sub$Healthy <- factor(Metadata_Df_sub$Healthy)
    
    table(Metadata_Df_sub$Country, useNA = 'ifany')
    
    RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_sub,
                           input_metadata = Metadata_Df_sub,
                           output = 'LRT_M13_AllBAL_Dis_Pls_All_ImpCov_output',
                           formula = '~ Healthy + (1 | Country) + Total_reads + observed',
                           normalization = 'TSS',
                           transform = 'LOG',
                           augment = TRUE,
                           standardize = TRUE,
                           max_significance = 0.1,
                           median_comparison_abundance = TRUE,
                           median_comparison_prevalence = FALSE,
                           max_pngs = 250,
                           cores = 2) 
  
    # Model 14 (M14) Effect on country by global dataset
    RT_fit_out <- maaslin3(input_data = Metaphlan_profiles_sub,
                           input_metadata = Metadata_Df_sub,
                           output = 'LRT_M14_AllBAL_Dis_Pls_All_ImpCov_output',
                           formula = '~ Healthy + Country + Total_reads + observed',
                           normalization = 'TSS',
                           transform = 'LOG',
                           augment = TRUE,
                           standardize = TRUE,
                           max_significance = 0.1,
                           median_comparison_abundance = TRUE,
                           median_comparison_prevalence = FALSE,
                           max_pngs = 250,
                           cores = 2) 

    Metadata_Df_DOI_OG <- Metadata_Df_DOI
  
    Metadata_Df_DOI <- Metadata_Df_DOI %>%
      filter(observed > 1)
      
  dim(Metadata_Df_DOI)
  
  table(Metadata_Df_DOI$BioProject, Metadata_Df_DOI$Disease)
  
  plot(Metadata_Df_DOI$BioProject, Metadata_Df_DOI$observed)
  