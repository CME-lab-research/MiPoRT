# Code to plot genome encoded functional potentials (KEGG module based) corresponding to KOs associated to disease in community profiles

# Load libraries
library(tidyverse)

# 1. Read input files

# Read filtered KEGG annotation file
KEGG_annotation <- 
    read_tsv("KEGG_annotations_for_MAGs_Maaslin-3_species.txt",
             name_repair = 'minimal')

    # Keep columns of interest and rm rest
    KEGG_annotation <- KEGG_annotation %>%
        select(all_of(c(
            "Species", "module", "module_name", "steps", "steps_present",
            "step_coverage", "ko_count", "kos_present"
        )))

    # sanity check
    glimpse(KEGG_annotation)
    summary(KEGG_annotation)
    
# Read Maaslin-3 file to get KO association details
LRT_associations <- read_tsv("LRT_KO_m1_associations_plus_m2_summary.txt",
                             name_repair = "minimal")

IRT_associations <- read_tsv("IRT_KO_m1_associations_plus_m2_summary.txt",
                             name_repair = "minimal")

    # Get KOs associated with Total reads and remove them
    KOs_to_rm <- unique(LRT_associations %>% 
        filter(Total_reads_flag == TRUE) %>%
        pull(feature))
    
        # rm such associations from the LRT table
        LRT_associations <- LRT_associations %>%
            filter(!feature %in% KOs_to_rm,
                   qval_individual < 0.05)%>%
            mutate(SampleType = "BAL")
    
    KOs_to_rm <- unique(IRT_associations %>% 
                                filter(Total_reads_flag == TRUE) %>%
                                pull(feature))
    
        # rm such associations from the IRT table
        IRT_associations <- IRT_associations %>%
            filter(!feature %in% KOs_to_rm,
                   qval_individual < 0.05) %>%
            mutate(SampleType = "Sputum")
    
    # Join/rbind LRT & IRT tables
        Columns_to_Keep <- c(
            "SampleType", "ModelID", "feature", "KO_Name", "value", "model",
            "coef", "qval_individual", "qval_joint")
        
        # rbind 2 tables
        Merged_association_df <- 
            rbind(
                LRT_associations[, Columns_to_Keep],
                IRT_associations[, Columns_to_Keep])
        
    # add groups and factors to associations table
    Merged_association_df <- Merged_association_df %>%
        mutate(model = 
                   factor(model,
                          levels = c('prevalence', 'abundance'),
                          labels = c('Prevalence', 'Abundance')
                          )
        )
    
    # sanity check
    table(Merged_association_df$value,
          Merged_association_df$model)

# 2. Format the KEGG annotation file
    glimpse(KEGG_annotation)
    summary(KEGG_annotation)
    
    glimpse(Merged_association_df)

    # sanity check
    # Check number of unique KOs to look for in genomes
    length(unique(Merged_association_df$feature))
        # 88
    
    # 2.1 Convert genome KEGG annotations into species-KO presence.
    
    # Split comma-separated KO annotations
    # Each kos_present cell contains multiple KOs separated by commas.
    # Convert these into one KO per row so they can be joined
    # directly to the MaAsLin3 association results.
    Genome_KO_presence <- KEGG_annotation %>%
        select(
            Species,
            kos_present
        ) %>%
        
        # Remove missing or empty KO annotations.
        filter(
            !is.na(kos_present),
            kos_present != ""
        ) %>%
        
        # Split the comma-separated KO string into separate rows.
        separate_longer_delim(
            kos_present,
            delim = ","
        ) %>%
        
        # Rename the column to match the biological entity.
        rename(
            KO = kos_present
        ) %>%
        
        # Remove any accidental white-space around KO identifiers.
        mutate(
            KO = stringr::str_trim(KO)
        ) %>%
        
        # A KO can occur in several KEGG modules for the same genome.
        # Keep each species-KO combination only once.
        distinct(
            Species,
            KO
        )
    
        # sanity check
        summary(Genome_KO_presence)
        
        # count KOs per genome
        Genome_KO_presence %>%
            count(
                Species,
                name = "Genome_KO_count"
            ) %>%
            arrange(
                desc(Genome_KO_count)
            )
    
    
    '
        # A tibble: 10 × 2
       species                     Genome_KO_count
       <chr>                                 <int>
     1 Veillonella atypica                     303
     2 Prevotella histicola                    272
     3 Leptotrichia wadei                      265
     4 Prevotella oris                         265
     5 Rothia mucilaginosa                     256
     6 Prevotella pallens                      246
     7 Prevotella jejuni                       244
     8 Porphyromonas endodontalis              233
     9 Peptostreptococcus stomatis             210
    '
    
# 3. Prepare significant MaAsLin3 KO associations
    
    # ---- 3.1 Keep unique disease-model-KO associations ----
    Association_KOs <- Merged_association_df %>%
        transmute(
            SampleType = SampleType,
            Disease = value,
            model = model,
            KO = feature,
            coef = coef
        ) %>%
        
        # Prevent duplicate association rows from causing
        # the same KO to be counted more than once.
        distinct(
            SampleType,
            Disease,
            model,
            KO,
            .keep_all = TRUE
        ) %>%
        mutate(
            Assoc_direction = if_else(coef > 0, 
                                      "Positive",
                                      "Negative"),
            Assoc_direction = factor(Assoc_direction),
            Disease = factor(Disease,
                             levels = c("Pneumonia", 
                                        "Covid-19",
                                        "Cystic fibrosis"))
        )
    
        # sanity check: Count unique KOs across the complete association table
        Association_KOs %>%
        summarise(
            Unique_KOs = n_distinct(KO)
        ) 
        
        # sanity check
        glimpse(Association_KOs)
        
        # Count KOs per model, Disease
        Association_KOs %>%
            count(
                SampleType,
                Disease,
                model,
                Assoc_direction,
                name = "Significant_KO_count"
            )

    '
# A tibble: 11 × 5
   SampleType Disease   model Assoc_direction Significant_KO_count
   <chr>      <fct>     <fct> <fct>                          <int>
 1 BAL        Pneumonia Prev… Negative                          31
 2 BAL        Pneumonia Prev… Positive                           2
 3 BAL        Pneumonia Abun… Negative                          25
 4 BAL        Pneumonia Abun… Positive                           6
 5 BAL        Covid-19  Prev… Negative                          11
 6 BAL        Covid-19  Prev… Positive                           9
 7 BAL        Covid-19  Abun… Negative                          14
 8 BAL        Covid-19  Abun… Positive                          14
 9 Sputum     Cystic f… Prev… Negative                           8
10 Sputum     Cystic f… Abun… Negative                           1
11 Sputum     Cystic f… Abun… Positive                          15
'
    
# 4. Intersect genome KOs with disease-associated KOs
    # Join genome-encoded KOs to MaAsLin3 associations
    # A row is retained only when:
    # 1. the KO is encoded by that species genome, and
    # 2. the KO is significantly associated with a disease/model.
    Genome_disease_KO_matches <- Genome_KO_presence %>%
        inner_join(
            Association_KOs,
            by = "KO",
            relationship = "many-to-many"
        )        

    # ---------- NOTE: Important KO-genomic KO mapping file ---------- 
    # save file
    write_tsv(Genome_disease_KO_matches,
              "Genome_KO_Corresponging_maaslin_KOs.txt")
    # ---------- Important KO-genomic KO mapping file ---------- 
    
    # ---- Count unique matched KOs ----
    Genome_disease_KO_counts <- Genome_disease_KO_matches %>%
        group_by(
            Species,
            SampleType,
            Disease,
            model,
            Assoc_direction
        ) %>%
        summarise(
                
            # Count each significant KO once per species.
            Disease_KO_count = n_distinct(KO),
                
            # Retain the actual KO identifiers for inspection.
            Associated_KOs = paste(
                sort(unique(KO)),
                collapse = ", "
            ),
                
            .groups = "drop"
        )        
        
    # sanity check
    Genome_disease_KO_counts
    
# 5. Add zero-count disease-model combinations. 
    # Define all observed disease-model comparisons
    # Use only combinations that actually occur in the
    # MaAsLin3 association results.
    Association_groups <- Association_KOs %>%
        distinct(
            SampleType,
            Disease,
            model,
            Assoc_direction
        )
    
    # Create all species-comparison combinations        
    # Every selected species should have one row for every
    # disease-model combination, including combinations with zero KOs.
    Complete_grid <- Association_groups %>%
        tidyr::crossing(
            Species = unique(Genome_KO_presence$Species)
        ) %>%
        relocate(
            Species,
            .before = SampleType
        )

# 6. Join observed counts onto the complete grid ----
    Disease_KO_summary <- Complete_grid %>%
        left_join(
            Genome_disease_KO_counts,
            by = c(
                "Species",
                "SampleType",
                "Disease",
                "model",
                "Assoc_direction"
            )
        ) %>%
        
        # Replace missing matches with zero.
        mutate(
            Disease_KO_count = tidyr::replace_na(
                Disease_KO_count,
                0L
            ),
            Associated_KOs = tidyr::replace_na(
                Associated_KOs,
                ""
            )
        )        
        

# 7. Add total significant-KOs count and percentages to the df
    # Count the total number of significant KOs identified
    # by MaAsLin3 for each disease-model combination.
    Association_KO_totals <- Association_KOs %>%
        group_by(
            SampleType,
            Disease,
            model,
            Assoc_direction
        ) %>%
        summarise(
            Total_significant_KOs = n_distinct(KO),
            .groups = "drop"
        )

    # Calculate genome representation percentage
    Disease_KO_summary <- Disease_KO_summary %>%
        left_join(
            Association_KO_totals,
            by = c(
                "SampleType",
                "Disease",
                "model",
                "Assoc_direction"
            )
        ) %>%
        mutate(
            
            # Percentage of disease-associated KOs represented
            # within the selected genome.
            Percent_encoded = 100 *
                Disease_KO_count /
                Total_significant_KOs
        )        

    # sanity check and sort
    Disease_KO_summary <- Disease_KO_summary %>%
        arrange(
            Species,
            Disease,
            model
        )
    
    # sanity check
    summary(Disease_KO_summary)
    
    # Add factors to relevant columns
    Disease_KO_summary <- Disease_KO_summary %>%
        mutate(
            Disease = factor(Disease, 
                             levels = c(
                                 "Pneumonia", "Covid-19", "Cystic fibrosis"
                             )),
            SampleType = factor(SampleType),
            model = factor(model)
        )
    
    summary(Disease_KO_summary)
    
    # Make a plotting df
    Disease_KO_summary_plot <- Disease_KO_summary %>%
        filter(
            Disease_KO_count > 0
        ) %>%
        droplevels() %>%
        
        # ---- 1.2 Set plotting order for disease groups ----
    
    mutate(
        Disease = factor(
            Disease,
            levels = c(
                "Pneumonia",
                "Covid-19",
                "Cystic fibrosis"
            )
        ),
        
        # This interaction defines how points are separated within each species.
        # It prevents disease/model/direction combinations from sitting exactly on top of each other.
        Plot_group = interaction(
            Disease,
            model,
            Assoc_direction,
            drop = TRUE
        )
    )
        
# 8. Plot disease-associated KO representation
    
    # Define shared plotting settings
    
    # Use the same dodge position for points and text
    disease_ko_dodge <- position_dodge(
            width = 0.45
        )
    
    # plot
    Disease_KO_plot <- ggplot(
            Disease_KO_summary_plot,
            aes(
                x = Species,
                y = Percent_encoded,
                shape = Disease,
                group = Plot_group
            )
        ) +
        
        # Show the % of disease-associated KOs encoded by each genome.
        geom_point(
            aes(
                color = model
            ),
            position = disease_ko_dodge,
            size = 5,
            alpha = 0.75,
            stroke = 0.25,
            na.rm = TRUE
        ) +
    
        # Print the number of encoded disease-associated KOs above each point.
        # This keeps both the normalized percentage and the raw KO count visible.
        geom_text(
            aes(
                label = Disease_KO_count
            ),
            position = disease_ko_dodge,
            vjust = -0.8,
            size = 6,
            show.legend = FALSE,
            na.rm = TRUE
        ) +
        
        # Separate positive and negative disease-associated KOs.
        # This makes association direction explicit instead of hiding it only in the grouping variable.
        #facet_wrap(
        #    ~ Assoc_direction,
        #    nrow = 1, 
        #    drop = TRUE
        #) +
        
        # Italicize species names on the x-axis.
        scale_x_discrete(
            labels = function(x) {
                parse(
                    text = paste0(
                        'italic("',
                        x,
                        '")'
                    )
                )
            }
        ) +
        
        # Show percentages as percentage values.
        scale_y_continuous(
            labels = scales::label_number(
                suffix = "%"
            ),
            expand = expansion(
                mult = c(
                    0.02,
                    0.15
                )
            )
        ) +
        
        # Add clean axis and legend labels.
        labs(
            x = NULL,
            y = "% disease-associated KOs\nrepresented in genome",
            title = "Genome-encoded representation of disease-associated community KOs",
            caption = "Point labels show the number of disease-associated HUMAnN KOs\nmatched to KEGG annotations encoded by the corresponding RESPICAT genome.",
            shape = "Disease associated",
            color = "Type of KO association"
        ) +
        
        # Increase size of legend icons
        guides(
            
            shape = 
                guide_legend(override.aes = list(size = 5)),
            color = 
                guide_legend(override.aes = list(size = 5)),
            ) +
    
        # Use a simple publication-style theme.
        theme_bw(
            base_size = 14
        ) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(),
            
            # Axis text-tick label settings
            axis.text.x = element_text(
                size = 14,
                angle = 45,
                hjust = 1,
                vjust = 1,
                color = "black"
            ),
            
            # Axis title settings
            axis.title = element_text(
                size = 16,
                face = "bold"
            ),
            
            # legend settings
            legend.text = element_text(size = 14,
                                       face = "bold",
                                       color = "grey20"),
            legend.title = element_text(
                size = 16,
                face = "bold"
            ),
            
            # legend position
            # Keep legends on the bottom
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vertical",
            legend.margin = margin(t = 2, r = 5, b = 2, l = 5, unit = "pt"),
            strip.text = element_text(
                face = "bold"
            )
        )

    # see it
    Disease_KO_plot    
    
    # save plot with less Taxa
    ggsave("./Plots/Genome-encoded_sp-KO_plot.png", 
           plot = Disease_KO_plot,
           width = 24,              # Reduce dimensions
           height = 24,
           units = "cm",
           dpi = 300,               # High resolution
           limitsize = FALSE
    )
    
    
    
    
    
    
    
    
    
    
    
        