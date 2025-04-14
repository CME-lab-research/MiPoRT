# Script to filter out singleton Taxa with a list of Taxa from PA table from the main Global Metaphlan4 species table. 

# load libs
library(tidyverse)

# read Data
FeatureTable_Df <- read.table("../MetaPhlan4_results/Metaphlan4_species_merged_filtered.tsv", header = T)

Filtered_list_Taxa <- row.names(read.table("../MetaPhlan4_results/MiPORT_filtered_PA_data.txt", header = T))
Filtered_list_Samples <- colnames(read.table("../MetaPhlan4_results/MiPORT_filtered_PA_data.txt", header = T))

# Generate filtered table with join
Filtered_FeatureTable_Df <- FeatureTable_Df %>%
    filter(taxonomy %in% Filtered_list_Taxa) %>%
    select(all_of(c("taxonomy", Filtered_list_Samples)))

# Write this table
write.table(Filtered_FeatureTable_Df, "../MetaPhlan4_results/Metaphlan4_species_filtered_Min2_samples.tsv", row.names = F, quote = F, sep = "\t")

# Write the Taxa and Samples names into a file as well?

Filtered_rowCol_Names <- list(
    Samples <- Filtered_list_Samples,
    Taxa <- Filtered_list_Taxa
)

    # collapse and write this list with paste and sapply
    
    tempDf <- data.frame(Names = sapply(Filtered_rowCol_Names, paste, collapse = ',')) 
    row.names(tempDf) <- c("Samples", "Species")
        
    write.table(tempDf, "../MetaPhlan4_results/Filter_Samples_Taxa.txt", sep = "\t", row.names = T)
