# Code to plot species prevalence in different healthy sample-types
#    - Add a layer of annotation/color if we we have recovered MAG for that species

# load libraries
library(tidyverse)

# 1. Read input files
    # 1.A: Read prevalence summary file
    prevalence_df <- read_tsv(
        "./Summary_PA_and_Abundance_Species_min10.txt",
        name_repair = "minimal")

    # subset columns
    Cols_to_Keep <- c("Species",
                      "Total_Prevalence-x",
                      "SamplePresence_count",
                      "Species_names",
                      colnames(prevalence_df) %>% 
                          str_subset("^(Prevalence_Healthy)")
                      )

    # select cols to keep
    prevalence_df <- prevalence_df %>%
        select(all_of(Cols_to_Keep)) %>%
        mutate(RT_prevalence = `Total_Prevalence-x`,
               Total_Prevalence.x = NULL) %>%
        relocate(RT_prevalence, .after = "SamplePresence_count")
    
    # format species names
    # Load the helper function
    source("C:/Users/o_shinde/Desktop/MiPORT_DA/Scripts/R_Utility_scripts/format_species_label.R")  # adjust path if needed
    
    # Format variable 'feature' and  create a new col 'Species names'
    prevalence_df$Species <- 
        format_species_label(
            prevalence_df$Species,
            italicize = FALSE)
    
    # 1.B: Read Maaslin-3 file
    association_df <- read_tsv("./Maaslin-3-associated_RESPICAT_intersected_Species.txt")
    
    # 1.3 Read RESPICAT taxonomy file
    RESPICAT_genome_taxonomy <- read_tsv(file = "./01_Select_species/RESPICAT_MAG_taxonomy_v1.0.tsv")
    
    # 1.4 Read RESPICAT genome features
    RESPICAT_genome_features <- 
        read_tsv(
            file = "./RESPICAT_MAG_quality_summary_v1.0.tsv",
            name_repair = "minimal")
    
    # Select genomes with species assignments
    RESPICAT_genome_taxonomy <- 
        RESPICAT_genome_taxonomy %>%
        filter(no_species_assignment == FALSE)
    
    # names of common species
    Species_to_Retain <- intersect(
        prevalence_df$Species,
        RESPICAT_genome_taxonomy$species
    )
    
        # Use this to filter the genome features df
        MAG_IDs_To_keep <- RESPICAT_genome_taxonomy %>%
            filter(species %in% Species_to_Retain) %>%
            distinct(species, .keep_all = TRUE) %>%
            pull(Public_MAG_ID)
    
    # 202 species from filtered M4 profiles have high-quality recovered genomes.
    
    # Check overlap between Maaslin-3 select species and prevalence select species
    intersect(Species_to_Retain, 
              unique(association_df$Species))
    
# 2. Filter prevalence df for species of interest (Top Maaslin-3, Top prevalence etc).
    
    # change this list value to update viz of select species
    Species_filter <- Species_to_Retain
    
    # OR
    # Filter-1: Select top 50 highly prevalent species in LRT >
    #   intersect with Maaslin-3 select species
    top_LRT_species <- prevalence_df %>%
        filter(Species %in% Species_to_Retain) %>%
        arrange(desc(Prevalence_Healthy_BAL)) %>%
        slice_head(n = 50) %>%
        pull(Species)
    
    # intersect with select Maaslin-3 species names
    Species_filter <- union(top_LRT_species,
                            unique(association_df$Species))

# 3. Subset data and prepare for plotting   
    # 3.1 Subset and melt the prevalence df
    plot_df <- prevalence_df %>%
        select(!all_of(
            c("SamplePresence_count",
              "RT_prevalence",
              "Total_Prevalence-x",
              "Prevalence_Healthy_Other",
              "Prevalence_Healthy_Other"
              ))) %>%
        filter(Species %in% Species_filter) %>%
        pivot_longer(
            cols = -c("Species", "Species_names"),
            names_to = "SampleType",
            values_to = "Species_Prevalence"
        ) %>%
        mutate(
            SampleType = str_remove(SampleType, "Prevalence_Healthy_")
        )
    # should ideally give you 56*7=392 rows
    
    # Add all factors to plot_df columns
    
    # init factors and labels
    SampleType_Factor <- c(
        "Nasal swab", 
        "Buccal mucosa", 
        "Oral swab", 
        "Saliva", 
        "Tongue dorsum", 
        "Sputum", 
        "BAL")
    
    # Add sampletype factor
    plot_df$SampleType <- factor(
        plot_df$SampleType,
        levels = SampleType_Factor
    )
    
        levels(plot_df$SampleType)
    
    # 3.2 # Initialize the plot with sample type on x-axis
        # ---- 3.1 Initialize prevalence heatmap ----
        GGplot_Obj <- ggplot(
            plot_df,
            aes(
                x = SampleType,
                y = Species,
                fill = Species_Prevalence
            )
        )
        
        
        # ---- 3.2 Draw prevalence tiles ----
        
        # Use a subtle border so neighbouring cells remain distinguishable
        # without making the grid visually dominant.
        GGplot_Obj <- GGplot_Obj +
            geom_tile(
                colour = "grey80",
                linewidth = 0.2
            ) +

        # ---- 3.3 Define prevalence colour scale ----
        # Keep prevalence on an explicit 0-100% sequential scale.
        scale_fill_gradientn(
                colours = c(
                    "white",
                    "#FEE8C8",
                    "#FDBB84",
                    "#FC8D59",
                    "#D7301F"
                ),
                limits = c(0, 100),
                breaks = c(0, 25, 50, 75, 100),
                labels = function(x) paste0(x, "%"),
                na.value = "grey90",
                name = "Prevalence"
            ) +
        
        
        # ---- 1.4 Control tile proportions ----
        # Keep the same compact tile geometry used previously.
        coord_fixed(
            ratio = 0.7
        ) +
        
        # ---- 1.5 Add concise publication-style labels ----
        # Shorten the title because the figure legend can provide
        # the full RESPICAT and healthy-sample context.
        labs(
            title = "Taxon prevalence across healthy respiratory sample types",
            x = NULL,
            y = "RESPICAT-represented taxa"
            ) +
        
        # ---- 1.6 Italicize species names ----
        scale_y_discrete(
            labels = function(x) {
                parse(
                    text = paste0('italic("', x, '")')
                )
            },
            expand = expansion(mult = c(0, 0))
        ) +

        # ---- 1.7 Remove unnecessary x-axis padding ----
        scale_x_discrete(
            expand = expansion(mult = c(0, 0))
        ) +

        # ---- 1.8 Apply publication-style theme ----
        theme_minimal(
            base_size = 11
        ) +
            theme(
            
            # Tile borders already provide visual structure.
            panel.grid = element_blank(),
                
            # Keep sample-type labels readable.
            axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1,
                size = 10,
                colour = "black"
            ),
                
            # Keep species labels clearly readable.
            axis.text.y = element_text(
                size = 10,
                colour = "black"
            ),
                
            # Make the y-axis title less visually dominant.
            axis.title.y = element_text(
                size = 11,
                face = "bold",
                margin = margin(r = 6)
            ),
                
            # Keep the title compact.
            plot.title = element_text(
                size = 13,
                face = "bold",
                hjust = 0
            ),
            
            plot.title.position = "plot",
            
            # Keep legend typography compact.
            legend.title = element_text(
                size = 10,
                face = "bold"
            ),
            
            legend.text = element_text(
                size = 9
            ),
                
            # Keep the legend below the heatmap.
            legend.position = "bottom",
            
            # Reduce outer whitespace so the genome-feature
            # block can sit closer to the heatmap.
            plot.margin = margin(
                t = 2,
                r = 0,
                b = 0,
                l = 2
            )
        ) +
        
        
        # ---- 1.9 Refine legend dimensions ----
        
        # Make the prevalence legend slightly wider and bring it
        # visually closer to the heatmap.
        guides(
            fill = guide_colourbar(
                title.position = "left",
                title.hjust = 0.5,
                barwidth = unit(5.5, "cm"),
                barheight = unit(0.35, "cm")
                )
        )
    
    # plot
    GGplot_Obj
    
    # 4. Prepare genome-feature annotation table
    RESPICAT_genome_features <- RESPICAT_genome_features %>%
        filter(Public_MAG_ID %in% MAG_IDs_To_keep)
    
    # Join taxonomy and feature dfs
    Merged_RESPICAT_feature_df <- 
        RESPICAT_genome_features %>%
        left_join(RESPICAT_genome_taxonomy,
                  by = join_by(Public_MAG_ID == Public_MAG_ID)) %>%
        relocate(species, .after = "Public_MAG_ID")
    
    glimpse(Merged_RESPICAT_feature_df)
    
    # Prepare plot-2 table
        # Keep only the species identifier and the two genome metrics
        # that will be displayed beside the prevalence heatmap.
        Genome_feature_plot_df <- 
            Merged_RESPICAT_feature_df %>%
            select(
                species,
                Size,
                N50
            ) %>%
            filter(
                species %in% Species_filter
            )

        # Reshape genome features into long format ----
        # Convert Size and N50 into a single feature column.
        # This gives one tile per species-feature combination.
        Genome_feature_plot_df <- 
                Genome_feature_plot_df %>%
                pivot_longer(
                    cols = c(Size, N50),
                    names_to = "Feature",
                    values_to = "Value"
                    ) %>%
                # Rename the feature columns exactly as they should appear
                mutate(
                    Feature = recode(
                        Feature,
                        "Size" = "Genome size",
                        "N50" = "N50"
                    )
                )
    
        # Create formatted text labels
        # Genome size is shown in megabases.
        # N50 is shown in kb unless it reaches at least 1 Mb.
        Genome_feature_plot_df <- Genome_feature_plot_df %>%
            mutate(
                Display_value = case_when(
                    
                    # Genome size is always shown in megabases.
                    Feature == "Genome size" ~ paste0(
                        round(Value / 1e6, 2),
                        " M"
                    ),
                    
                    # N50 values >= 1,000,000 bp are shown in megabases.
                    Feature == "N50" & Value >= 1e6 ~ paste0(
                        round(Value / 1e6, 2),
                        " M"
                    ),
                    
                    # Smaller N50 values are shown in kilobases.
                    Feature == "N50" ~ paste0(
                        round(Value / 1e3, 1),
                        " K"
                    )
                )
            )
    
        # Keep genome size first and N50 second.
        Genome_feature_plot_df <- Genome_feature_plot_df %>%
                mutate(
                    Feature = factor(
                        Feature,
                        levels = c(
                            "Genome size",
                            "N50"
                        )
                    )
                )
            
        # Assign Genome size and N50 to positions with a gap between them.
        Genome_feature_plot_df <- Genome_feature_plot_df %>%
            mutate(
                Feature_position = case_when(
                    Feature == "Genome size" ~ 1,
                    Feature == "N50" ~ 2.15
                )
            )
            

        # Build genome-feature annotation plot
        # ---- Initialize the plot ----
        Genome_feature_plot <- ggplot(
            Genome_feature_plot_df,
            aes(
                x = Feature_position,
                y = species
            )
        )
            
            
        # ---- Draw neutral tiles ----
            
        # Use a light neutral fill because the actual quantitative
        # information is communicated directly by the text labels.
        Genome_feature_plot <- Genome_feature_plot +
                geom_tile(
                    width = 0.85,
                    height = 0.95,
                    fill = "grey96",
                    colour = "grey82",
                    linewidth = 0.35
                ) + 
                coord_fixed(ratio = 0.7) +
            
        # ---- Add formatted genome values ----
        # Print the genome-size and N50 values directly inside each tile.
        geom_text(
                aes(
                    label = Display_value
                    ),
                size = 3.0,
                colour = "black"
            ) +
            
            # Show labels only at the positions containing tiles.
            scale_x_continuous(
                breaks = c(1, 2.15),
                labels = c(
                    "Genome size",
                    "N50"
                ),
                expand = expansion(
                    mult = c(0.01, 0.01)
                )
            ) +
            
            # Italicize labels
            scale_y_discrete(
                labels = function(x) {
                    parse(
                        text = paste0(
                            'italic("',
                            x,
                            '")'
                        )
                    )
                },
                drop = FALSE,
                expand = expansion(mult = c(0, 0))
            )  +
        
            
            # apply publication-style formatting
            # ---- add labels ----
            labs(
                x = NULL,
                subtitle = "Genome features",
                y = NULL
            ) +
            
            
            # ---- Apply a clean theme ----
            theme_minimal(
                base_size = 10
                ) +
            theme(
                # Remove grid lines because tile boundaries
                # already define the table structure.
                panel.grid = element_blank(),
                
                # Keep feature labels readable.
                axis.text.x = element_text(
                    size = 9.5,
                    angle = 90,
                    hjust = 1,
                    vjust = 0.5,
                    colour = "black",
                    face = "bold",
                    margin = margin(t = 4)
                    ),
                
                # Keep species labels visible for now.
                axis.text.y = element_text(
                    size = 10,
                    colour = "black"
                    ),
                
                # Remove axis ticks.
                axis.ticks = element_blank(),
                    
                # Reduce unnecessary outer whitespace.
                plot.margin = margin(
                        t = 2,
                        r = 2,
                        b = 0,
                        l = 0
                    )
                
                )
            
        # Display genome feature plot
        Genome_feature_plot 
        
        # save png plot
        ggsave("./Plots/RESPICAT_genomic_features.png", 
               plot = Genome_feature_plot,
               width = 18,              # Reduce dimensions
               height = 48,
               units = "cm",
               dpi = 600,               # High resolution
               limitsize = FALSE
        )
            
        # save svg plot
        ggsave("./Plots/RESPICAT_genomic_features.svg",
               plot = Genome_feature_plot,
               width = 18,              # Reduce dimensions
               height = 48,
               units = "cm",
               limitsize = FALSE
            )
    
    
    ### combine plots
    
    # Keep species labels only on the prevalence heatmap.
    Genome_feature_plot_clean <- Genome_feature_plot +
        theme(
            axis.text.y = element_blank()
        )
            
            
    # ---- 1.2 Combine panels with controlled widths ----
    
    # Give the prevalence matrix substantially more horizontal space
    # than the two-column genomic annotation table.
    Combined_plot <- GGplot_Obj +
        Genome_feature_plot_clean +
        patchwork::plot_layout(
            widths = c(7, 3)
        )
            
    Combined_plot
    
    # save png plot
    ggsave("./Plots/RESPICAT_heatmap_integration_plot.png", 
           plot = Combined_plot,
           width = 25,              # Reduce dimensions
           height = 45,
           units = "cm",
           dpi = 600,               # High resolution
           limitsize = FALSE
    )
    
    ggsave("./Plots/RESPICAT_heatmap_integration_plot.svg", 
           plot = Combined_plot,
           width = 25,              # Reduce dimensions
           height = 45,
           units = "cm",
           limitsize = FALSE
    )
    
    
    



