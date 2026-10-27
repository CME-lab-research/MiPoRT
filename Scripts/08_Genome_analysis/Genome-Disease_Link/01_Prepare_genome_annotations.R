# Code to filter the table with KEGG output for all genomes (long format)
# Input file: RESPICAT_MAGs_dRep_KEGG_output_long_df.txt


# load library
library(tidyverse)

# 1. Read input files

# Read KEGG annotation table
All_KEGG_annotations <- read_tsv("RESPICAT_MAGs_dRep_KEGG_output_long_df.txt",
                                 name_repair = "minimal")

# Read species of interest file
Species_list <- read_tsv("Maaslin-3_Based_List_of_Species.txt", 
                         name_repair = "minimal",
                         col_names = FALSE
                        )
    # sort list
    Species_list <- Species_list %>% 
        arrange(desc(X1))
    # sanity check
    Species_list
    
# Read RESPICAT genome mapping file
RESPICAT_genome_mapping_file <- 
read_tsv("RESPICAT_files/RESPICAT_MAG_taxonomy_v1.0.tsv",
         name_repair = "minimal")

    # filter table with missing species annotations
    RESPICAT_genome_mapping_file <- RESPICAT_genome_mapping_file %>%
        filter(no_species_assignment == FALSE) %>%
        select(-no_species_assignment)
    
# Read the 2nd mapping file
RESPICAT_genome_ID_mapping_file <- 
    read_tsv("RESPICAT_files/RESPICAT_MAG_id_mapping_v1.0.tsv",
             name_repair = "minimal")
    
    # filter based on above mapping file
    RESPICAT_genome_ID_mapping_file <- RESPICAT_genome_ID_mapping_file %>%
        filter(Public_MAG_ID %in% 
                   RESPICAT_genome_mapping_file$Public_MAG_ID) %>%
        select(all_of(c("Public_MAG_ID", "Unique_MAG_ID")))
    
# 2. Select for genomes in Species list
    Species_df <- Species_list %>% 
        left_join(
            RESPICAT_genome_mapping_file,
            by = join_by("X1" == "species"),
            multiple = "any"
            ) %>%
        mutate(Species = X1,
               X1 = NULL) %>%
        left_join(
            RESPICAT_genome_ID_mapping_file,
            by = join_by("Public_MAG_ID" == "Public_MAG_ID")
        )

# 3. Filter KEGG annotation file
    Filtered_KEGG_annotations <- All_KEGG_annotations %>%
        filter(Unique_MAG_ID %in% Species_df$Unique_MAG_ID) %>%
        left_join(Species_df %>% 
                      select(all_of(c("Unique_MAG_ID", "Species"))),
                  by = join_by("Unique_MAG_ID" == "Unique_MAG_ID")
                  ) %>%
        relocate(Unique_MAG_ID, 
                 .after = last_col()) %>%
        relocate(Species)
    
    # sanity check to count number of genomes/species
    table(Filtered_KEGG_annotations$Species)
    
# 4. save filtered annotation file
    write.table(Filtered_KEGG_annotations,
                "KEGG_annotations_for_MAGs_Maaslin-3_species.txt",
                sep = "\t",
                row.names = FALSE)
    
    
    
    