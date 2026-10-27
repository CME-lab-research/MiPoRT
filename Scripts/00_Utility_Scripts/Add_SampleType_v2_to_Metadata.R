# Add sampleType v2 to metadata

library(tidyverse)

# Read metadata file
Metadata_Df <- read.table("./MiPORT_Metadata_downstream_filtered_M4_passed.tsv", header = T, sep = "\t")

# 1. Add factors Sample type
SamplingSite_Factor <- c("Nasal_Swab", "Nasopharyngeal_Aspirate", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Palatine_Tonsils","Throat", "Sputum", "BAL")

SamplingSite_Labels <- c("Nasal swab", "Nasopharyngeal aspirate", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Palatine tonsils","Throat", "Sputum", "BAL")

# Add factors
Metadata_Df$SampleType <- factor(Metadata_Df$SampleType,
                                 levels = SamplingSite_Factor,
                                 labels = SamplingSite_Labels)

table(Metadata_Df$SampleType, useNA = 'ifany')

# Count frequency of each SampleType
Top_sampletypes <- Metadata_Df %>%
    count(SampleType, sort = TRUE) %>%   # Count and sort
    slice_max(n, n = 8) %>%              # Get top 8 most frequent
    pull(SampleType)    %>%                 # Extract SampleType names
    droplevels()

# Add a new SampleType col which only includes top 8 sampletypes
# Re-assign remaining sample types as "Other"
Metadata_Df_mod <- Metadata_Df %>%
    mutate(SampleTypev2 = ifelse(as.character(SampleType) %in% Top_sampletypes, as.character(SampleType), "Other"))

Metadata_Df_mod <- Metadata_Df_mod %>%
    select(c("SampleID", "BioProject", "SampleType", "SampleTypev2", "Disease", "Healthy", "Age", "AgeGroup", "Antibiotics", "RT_category", "PatientID", "Time_pt", "Gender", "Smoking", "Enrichment", "Country", "Comments", "ProjectID", "Before_QC_R1", "Before_QC_R2", "After_decontam_R1", "After_decontam_R2", "After_QC_R1", "After_QC_R2", "After_decontam_R1_percent",,"After_decontam_R2_percent", "QC_Status_R1", "QC_Status_R2", "dominance_gini", "observed", "diversity_shannon", "diversity_gini_simpson"))
    
# sanity check
    table(Metadata_Df_mod$SampleTypev2, useNA = 'ifany')

# Add factors
Metadata_Df_mod$SampleTypev2 <- factor(Metadata_Df_mod$SampleTypev2, levels = c(levels(Top_sampletypes), "Other"))  

# save the table
write.table(Metadata_Df_mod, 
            "./MiPORT_Metadata_downstream_filtered_M4_passed.tsv", 
            row.names = F, sep = "\t", quote = F)

















