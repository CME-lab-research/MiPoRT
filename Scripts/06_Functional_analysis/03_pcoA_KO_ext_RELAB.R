# Script to calculate beta diversity and ordinate the funtional profiles

# Load env
library(tidyverse)
library(vegan)
library(ape)

# Step-1: Read merged HUMANN abundance table and Metadata file
    KO_ext_relab <- read.table(
        "../03_BC_Humann4_merged_KO_extended_relab.tsv",
        sep = '\t',
        header = T,
        check.names = FALSE)
    
    # read metadata
    Metadata_df <- read.table(
        "../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv",
        sep = '\t',
        header = T,
        check.names = FALSE)

    # sanity check
    Common_samples <- 
        intersect(
            Metadata_df$SampleID, 
            colnames(KO_ext_relab)[-1])
    
    setdiff(Metadata_df$SampleID, colnames(KO_ext_relab))

    # subset Metadata df 
    Metadata_df <- Metadata_df %>%
        filter(SampleID %in% Common_samples) %>%
        select(all_of(c("SampleID", "BioProject", "SampleType", 
                        "SampleTypev2", "RT_category")))
    
    # Add RT factors
    table(Metadata_df$RT_category, useNA = 'ifany')
    
    RT_Cat_Factor <- c("URT", "IRT", "LRT")
    RT_Cat_Labels <- c("Upper RT (N=1630)", "Intermediate RT (N=719)", "Lower RT (N= 578)")
    
    # add factors
    Metadata_df$RT_category <- 
        factor(Metadata_df$RT_category,
               levels = RT_Cat_Factor,
               labels = RT_Cat_Labels
        )
    
    # sanity check
    table(Metadata_df$RT_category, useNA = 'ifany')
    
    # Add sampletype factor
    table(Metadata_df$SampleTypev2)
    
    # Factors
    SampleType_Factor <- c("Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Sputum", "BAL", "Other")
    
    # Update labels with sample size
    SampleType_Labels <- c(
        "Nasal swab (n=27)",
        "Buccal mucosa (n=1027)",
        "Oral swab (n=246)",
        "Saliva (n=315)",
        "Tongue dorsum (n=418)",
        "Supraglottal (n=59)",
        "Sputum (n=193)",
        "BAL (n=578)",
        "Other (n=64)"
    )
    
    # add factors
    Metadata_df$SampleTypev2 <- 
        factor(Metadata_df$SampleTypev2,
               levels = SampleType_Factor,
               labels = SampleType_Labels
        )

    # sanity check
    table(Metadata_df$SampleTypev2)
    
    # subset KO_ext_relab
    KO_ext_relab <- KO_ext_relab %>%
        select(all_of(c("KO_ext_relab", Common_samples)))

# Step-2: Prepare ordination data
    
    # 2.1: Transpose (vegan expects samples as rows)
    KO_ext_relab_t <- KO_ext_relab %>% 
        column_to_rownames("KO_ext_relab") %>%
        t()

    # 2.2 Prepare weights for weighted CMD
    
    # Calculate weights according to sample-type counts
    # sanity check
    table(Metadata_df$SampleType, useNA = 'ifany')
    
    # Total sample size per SampleType
    sampletype_counts <- table(Metadata_df$SampleType, 
                               useNA = 'ifany')
    
    # Calculate weights based on the reciprocal of sample size
    weights <- case_when(
        Metadata_df$SampleType == "BAL" ~ 1 / sampletype_counts["BAL"],
        Metadata_df$SampleType == "Buccal_mucosa" ~ 1 / sampletype_counts["Buccal_mucosa"],
        Metadata_df$SampleType == "Nasal_Swab" ~ 1 / sampletype_counts["Nasal_Swab"],
        Metadata_df$SampleType == "Nasopharyngeal_Aspirate" ~ 1 / sampletype_counts["Nasopharyngeal_Aspirate"],
        Metadata_df$SampleType == "Oral_swab" ~ 1 / sampletype_counts["Oral_swab"],
        Metadata_df$SampleType == "Palatine_Tonsils" ~ 1 / sampletype_counts["Palatine_Tonsils"],
        Metadata_df$SampleType == "Saliva" ~ 1 / sampletype_counts["Saliva"],
        Metadata_df$SampleType == "Sputum" ~ 1 / sampletype_counts["Sputum"],
        Metadata_df$SampleType == "Supraglottal" ~ 1 / sampletype_counts["Supraglottal"],
        Metadata_df$SampleType == "Throat" ~ 1 / sampletype_counts["Throat"],
        Metadata_df$SampleType == "Tongue_dorsum" ~ 1 / sampletype_counts["Tongue_dorsum"],
        TRUE ~ NA_real_  # Catch any unexpected values
    )
    
# Step-3: Calculate BC distance
    BC_dist <- vegdist(KO_ext_relab_t, 
                       method = "bray")
    
# Step-4: PCA
    nComponents <- 3
    pcoa_result_weighted <- wcmdscale(BC_dist,
                             k = nComponents,
                             eig = TRUE, 
                             w = weights # add sample type weights here
                             )
    
    # Also calculate unweighted pcoa
    pcoa_result <- cmdscale(BC_dist,
                             k = nComponents,
                             eig = TRUE)

    # Calculate the variance explained by each PC & Convert to percentages
    variance_explained_weighted <- 
        pcoa_result_weighted$eig / sum(pcoa_result_weighted$eig) * 100  
    
    # unweighted var
    variance_explained <- 
        pcoa_result$eig / sum(pcoa_result$eig) * 100  
    
    # Add variance explained to PC names 
    PC_names_w <- c(
        paste0(
            "Component ",
            1:nComponents,
            " (",
            round(variance_explained_weighted[1:nComponents], 2), 
            "%)")
    )
    
    # unweighted
    PC_names <- c(
        paste0(
            "Component ",
            1:nComponents,
            " (",
            round(variance_explained[1:nComponents], 2), 
            "%)")
    )

# Extract the coordinates for the PCs
    pcoa_coords <- as.data.frame(pcoa_result$points)
    colnames(pcoa_coords) <- PC_names
    
    # for weighted
    pcoa_coords_w <- as.data.frame(pcoa_result_weighted$points)
    colnames(pcoa_coords_w) <- PC_names_w

# Plot 1: Un-weighted PCoA
    # left join data matrix and sample IDs
    Merged_df <- pcoa_coords %>% 
    rownames_to_column(var="SampleID") %>% 
    left_join(Metadata_df, 
              by = join_by("SampleID" == "SampleID"))

    # glimpse factors
    table(Merged_df$SampleTypev2, 
          Merged_df$RT_category,
          useNA = 'ifany')
    
    # plot obj
    ggplotObj <- ggplot(Merged_df, 
                    aes(x = `Component 1 (28.42%)`, 
                        y= `Component 2 (10.37%)`, 
                        colour = RT_category))

    savePlot <- ggplotObj + 
    geom_point(size = 1.5, 
               position = position_jitter(width = 0.02)) + 
    labs(title = 'PCoA based on KO abundances',
         subtitle = "(Bray-Curtis distance)",
         colour = "RT category") + 
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "top",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
            ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
               )
    
    savePlot

    # save png
    ggsave("./Plots/PCoA_filtered_BC_KO_ext_RT_Category.png",
           plot = savePlot,
           height = 14,
           width = 26,          
           units = "cm",
           dpi = 600,               # High resolution
           limitsize = FALSE)
    
    # save svg
    ggsave("./Plots/PCoA_filtered_BC_KO_ext_RT_Category.svg",
           plot = savePlot,
           height = 14,
           width = 26,
           units = "cm")
  
# Plot 2: Weighted PCoA
    # left join data matrix and sample IDs
    Merged_df_weighted <- pcoa_coords_w %>% 
        rownames_to_column(var="SampleID") %>% 
        left_join(Metadata_df, 
                  by = join_by("SampleID" == "SampleID"))
    
    
    # glimpse factors
    table(Merged_df_weighted$SampleTypev2, 
          Merged_df_weighted$RT_category,
          useNA = 'ifany')
    
    # plot obj
    ggplotObj <- ggplot(Merged_df_weighted, 
                        aes(x = `Component 1 (19.9%)`, 
                            y= `Component 2 (11.65%)`, 
                            colour = RT_category))
    
    savePlot <- ggplotObj + 
        geom_point(size = 1.5, 
                   position = position_jitter(width = 0.02)) + 
        labs(title = 'PCoA based on KO abundances',
             subtitle = "(Bray-Curtis distance)",
             colour = "RT category") + 
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "top",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot
    
    # save png
    ggsave("./Plots/w_PCoA_filtered_BC_KO_ext_RT_Category.png",
           plot = savePlot,
           height = 14,
           width = 26,          
           units = "cm",
           dpi = 600,               # High resolution
           limitsize = FALSE)
    
    # save svg
    ggsave("./Plots/w_PCoA_filtered_BC_KO_ext_RT_Category.svg",
           plot = savePlot,
           height = 14,
           width = 26,
           units = "cm")
    
    # Same plot with ellipses
    # plot obj
    ggplotObj <- ggplot(Merged_df_weighted, 
                        aes(x = `Component 1 (19.9%)`, 
                            y= `Component 2 (11.65%)`, 
                            colour = RT_category))
    
    savePlot <- ggplotObj + 
        geom_point(size = 1.5, 
                   position = position_jitter(width = 0.02)) + 
        stat_ellipse(aes(fill = RT_category), 
                     geom = "polygon", 
                     alpha = 0.09, 
                     level = 0.95) +  # Confidence ellipses
        labs(title = 'PCoA based on KO abundances',
             subtitle = "(Bray-Curtis distance)",
             colour = "RT category") + 
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "top",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(fill = "none",
               colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot
    
    # save png
    ggsave("./Plots/w_ePCoA_filtered_BC_KO_ext_RT_Category.png",
           plot = savePlot,
           height = 16,
           width = 26,          
           units = "cm",
           dpi = 600,               # High resolution
           limitsize = FALSE)
    
    # save svg
    ggsave("./Plots/w_ePCoA_filtered_BC_KO_ext_RT_Category.svg",
           plot = savePlot,
           height = 14,
           width = 26,
           units = "cm")
    
    
    
# Plot 3: Weighted PCoA with sampletype colors
    
    # Define custom colors for top 8 + "Other"
    SampleType_custom_colors <- c(
        "Nasal swab (n=27)" = "#DAA520", 
        "Buccal mucosa (n=1027)" = "#F781BF",
        "Oral swab (n=246)" = "#FF7F00",
        "Saliva (n=315)" = "#A65628",
        "Tongue dorsum (n=418)" = "#4DAF4A",
        "Supraglottal (n=59)" = "#984EA3",
        "Sputum (n=193)" = "#E41A1C", 
        "BAL (n=578)" = "#377EB8",
        "Other (n=64)" = "#999999"  # Gray for "Other"
    )
    
    # glimpse factors
    table(Merged_df_weighted$SampleTypev2, 
          Merged_df_weighted$RT_category,
          useNA = 'ifany')
    
    # plot obj
    ggplotObj <- ggplot(Merged_df_weighted, 
                        aes(x = `Component 1 (19.9%)`, 
                            y= `Component 2 (11.65%)`, 
                            colour = SampleTypev2))
    
    # plot
    savePlot <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) + 
        scale_colour_manual(values = SampleType_custom_colors) +
        labs(title = 'PCoA based on KO abundances',
             subtitle = "(Bray-Curtis distance)",
             colour = "Sample type (N)") + 
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "top",
            #legend.position = "bottom",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot
    
    # save png
    ggsave("./Plots/w_PCoA_filtered_BC_KO_ext_SampleType.png",
           plot = savePlot,
           height = 16,
           width = 26,          
           units = "cm",
           dpi = 600,               # High resolution
           limitsize = FALSE)
    
    # save svg
    ggsave("./Plots/w_PCoA_filtered_BC_KO_ext_SampleType.svg",
           plot = savePlot,
           height = 16,
           width = 26,
           units = "cm")
   
    # same plot with ellipses
    # plot obj
    ggplotObj <- ggplot(Merged_df_weighted, 
                        aes(x = `Component 1 (19.9%)`, 
                            y= `Component 2 (11.65%)`, 
                            colour = SampleTypev2))
    
    # plot
    savePlot <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) + 
        scale_colour_manual(values = SampleType_custom_colors) +
        labs(title = 'PCoA based on KO abundances',
             subtitle = "(Bray-Curtis distance)",
             colour = "Sample type (N)") + 
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "top",
            #legend.position = "bottom",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(fill = 
                   guide_legend(override.aes = list(size = 5)),
               colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot
    
    # save png
    ggsave("./Plots/w_ePCoA_filtered_BC_KO_ext_SampleType.png",
           plot = savePlot,
           height = 16,
           width = 26,          
           units = "cm",
           dpi = 600,   # High resolution
           limitsize = FALSE)
    
    # save svg
    ggsave("./Plots/w_ePCoA_filtered_BC_KO_ext_SampleType.svg",
           plot = savePlot,
           height = 16,
           width = 26,
           units = "cm")
    
    
############################################################
####    3 panel plot for RT category color             #####
############################################################
    
    # Plot 2.1
    # plot PC1 vs PC2
    ggplotObj <- Merged_df_weighted %>% 
        ggplot(
            aes(x=`Component 1 (19.9%)`, 
                y=`Component 2 (11.65%)`, 
                color= RT_category
            )
        )
    
    # save plot
    savePlot <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        labs(colour = "Respiratory category (N)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    # plot PC2 vs PC3
    ggplotObj <- Merged_df_weighted %>%
        ggplot(aes(x=`Component 3 (7.69%)`, 
                   y=`Component 2 (11.65%)`, 
                   color=RT_category #, shape = outlier
        )
        )
    # save plot
    savePlot_1 <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        labs(colour = "Respiratory category (N)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot_1
    
    # plot PC1 vs PC3
    ggplotObj <- Merged_df_weighted %>% 
        ggplot(aes(x= `Component 3 (7.69%)`, 
                   y= `Component 1 (19.9%)`, 
                   color= RT_category
        )
        )
    
    savePlot_2 <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        labs(colour = "Respiratory category (N)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot_2
    
    # To save multiple plots
    # Load library 
    library(patchwork)
    
    # Combine the plots into a panel (e.g., 1 row with 3 columns)
    combined_plot <- 
        savePlot + 
        savePlot_1 + 
        savePlot_2 +
        plot_layout(
            ncol = 3,
            guides = "collect"
        ) +
        plot_annotation(
            title = "PCoA based on KO abundances",
            subtitle = "Bray-Curtis distance",
            tag_levels = "A"
        ) &
        theme(
            legend.position = "bottom",
            plot.title = 
                element_text(size = 30, 
                             face = "bold", 
                             hjust = 0.5),
            plot.subtitle = element_text(size = 26, 
                                         hjust = 0.5)
        )
    
    combined_plot
    
    # Save combined plot
    ggsave("Plots/PCoA_filtered_BC_KO_ext_combined_panel_RT_Category.png",
           plot = combined_plot,
           height = 25, 
           width = 56, 
           units = "cm", 
           dpi = 600)
    
    # save svg
    ggsave("./Plots/PCoA_filtered_BC_KO_ext_combined_panel_RT_Category.svg",
           plot = combined_plot,
           height = 25,
           width = 56,
           units = "cm")
    
############################################################
#### Repeat the 3 panel plot for RT category color now #####
############################################################
    
    # Plot 2.1
    # plot PC1 vs PC2
    ggplotObj <- Merged_df_weighted %>% # filter(outlier == "FALSE") %>% 
        ggplot(
            aes(x=`Component 1 (19.9%)`, 
                   y=`Component 2 (11.65%)`, 
                   color= RT_category
                )
        )
    
    # save plot
    savePlot <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        labs(colour = "Respiratory category (N)") +
        stat_ellipse(aes(fill = RT_category), 
                     geom = "polygon", 
                     alpha = 0.2, 
                     level = 0.95) +  # Confidence ellipses
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(fill = 'none',
               colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    # plot PC2 vs PC3
    ggplotObj <- Merged_df_weighted %>%
        ggplot(aes(x=`Component 3 (7.69%)`, 
                   y=`Component 2 (11.65%)`, 
                   color=RT_category #, shape = outlier
        )
        )
    # save plot
    savePlot_1 <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        labs(colour = "Respiratory category (N)") +
        stat_ellipse(aes(fill = RT_category), 
                     geom = "polygon", 
                     alpha = 0.2, 
                     level = 0.95) +  # Confidence ellipses
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(fill = 'none',
               colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot_1
    
    # plot PC1 vs PC3
    ggplotObj <- Merged_df_weighted %>% 
        ggplot(aes(x= `Component 3 (7.69%)`, 
                   y= `Component 1 (19.9%)`, 
                   color= RT_category
        )
        )
    
    savePlot_2 <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        labs(colour = "Respiratory category (N)") +
        stat_ellipse(aes(fill = RT_category), 
                     geom = "polygon", 
                     alpha = 0.2, 
                     level = 0.95) +  # Confidence ellipses
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(fill = 'none',
               colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot_2
    
# To save multiple plots
    # Load library 
    library(patchwork)
    
    # Combine the plots into a panel (e.g., 1 row with 3 columns)
    combined_plot <- 
        savePlot + 
        savePlot_1 + 
        savePlot_2 +
        plot_layout(
            ncol = 3,
            guides = "collect"
        ) +
        plot_annotation(
            title = "PCoA based on KO abundances",
            subtitle = "Bray-Curtis distance",
            tag_levels = "A"
        ) &
        theme(
            legend.position = "bottom",
            plot.title = 
                element_text(size = 30, 
                             face = "bold", 
                             hjust = 0.5),
            plot.subtitle = element_text(size = 26, 
                                         hjust = 0.5)
        )
    
    combined_plot
    
    # Save combined plot
    ggsave("Plots/PCoA_with Ellipses_filtered_BC_KO_ext_combined_panel_RT_Category.png",
           plot = combined_plot,
           height = 25, 
           width = 56, 
           units = "cm", 
           dpi = 600)
    
    # save svg
    ggsave("./Plots/PCoA_with Ellipses_filtered_BC_KO_ext_combined_panel_RT_Category.svg",
           plot = combined_plot,
           height = 25,
           width = 56,
           units = "cm")
    
    
############################################################
####    3 panel plot for sample type color             #####
############################################################
    
    # Plot 2.1
    # plot PC1 vs PC2
    ggplotObj <- Merged_df_weighted %>% 
        ggplot(
            aes(x=`Component 1 (19.9%)`, 
                y=`Component 2 (11.65%)`, 
                color= SampleTypev2
            )
        )
    
    # save plot
    savePlot <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        scale_colour_manual(values = SampleType_custom_colors) +
        labs(colour = "Sample type (N)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    # plot PC2 vs PC3
    ggplotObj <- Merged_df_weighted %>%
        ggplot(aes(x=`Component 3 (7.69%)`, 
                   y=`Component 2 (11.65%)`, 
                   color= SampleTypev2
        )
        )
    # save plot
    savePlot_1 <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        scale_colour_manual(values = SampleType_custom_colors) +
        labs(colour = "Sample type (N)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot_1
    
    # plot PC1 vs PC3
    ggplotObj <- Merged_df_weighted %>% 
        ggplot(aes(x= `Component 3 (7.69%)`, 
                   y= `Component 1 (19.9%)`, 
                   color= SampleTypev2
        )
        )
    
    savePlot_2 <- ggplotObj + 
        geom_point(size = 1.5,
                   position = position_jitter(width = 0.02)) +
        scale_colour_manual(values = SampleType_custom_colors) +
        labs(colour = "Sample type (N)") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 22), 
            axis.title = element_text(size = 24, face = "bold"),
            title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.title.position = "left",
            legend.key.size = unit(1, 'cm'),
            legend.spacing.y = unit(0.6, 'cm'),
            plot.margin = 
                margin(t = 20, 
                       r = 20, 
                       b = 20, 
                       l = 20, 
                       unit = "pt"),
            plot.subtitle = element_text(size = 22)
        ) +
        guides(colour = 
                   guide_legend(override.aes = list(size = 5),
                                byrow = FALSE)
        )
    
    savePlot_2
    
    # To save multiple plots
    # Load library 
    library(patchwork)
    
    # Combine the plots into a panel (e.g., 1 row with 3 columns)
    combined_plot_sampleType <- 
        savePlot + 
        savePlot_1 + 
        savePlot_2 +
        plot_layout(
            ncol = 3,
            guides = "collect"
        ) +
        plot_annotation(
            title = "PCoA based on KO abundances",
            subtitle = "Bray-Curtis distance",
            tag_levels = "A"
        ) &
        theme(
            legend.position = "bottom",
            plot.title = 
                element_text(size = 30, 
                             face = "bold", 
                             hjust = 0.5),
            plot.subtitle = element_text(size = 26, 
                                         hjust = 0.5)
        )
    
    combined_plot_sampleType
    
    # Save combined plot
    ggsave("Plots/PCoA_filtered_BC_KO_ext_combined_panel_SampleType.png",
           plot = combined_plot_sampleType,
           height = 25, 
           width = 56, 
           units = "cm", 
           dpi = 600)
    
    # save svg
    ggsave("./Plots/PCoA_filtered_BC_KO_ext_combined_panel_SampleType.svg",
           plot = combined_plot_sampleType,
           height = 25,
           width = 56,
           units = "cm")

    
# save 3D plot for sampletype
    # Load required library
    library(htmlwidgets)
    library(plotly)
     
    # Create the sample type plot
    ord_3d_plot <- plot_ly(
        data = Merged_df_weighted,
        x = ~ `Component 1 (19.9%)`,
        y = ~ `Component 2 (11.65%)`,
        z = ~ `Component 3 (7.69%)`,
        color = ~ SampleTypev2,
        colors = SampleType_custom_colors,
        text = ~ paste("SampleID:", SampleID,
                      "<br>RT Category:", RT_category,
                      "<br>Sample Type:", SampleTypev2),
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 4, alpha = 0.5)
    ) %>%
        layout(
            title = list(
                text = "Bray-Curtis PCoA based on KO abundances",
                font = list(size = 22, 
                            family = "Arial", 
                            color = "black")
            ),
            margin = list(t = 80), 
            scene = list(
                xaxis = list(
                    title = list(text = "<br><br>Component 1 (19.9%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                ),
                yaxis = list(
                    title = list(text = "<br><br>Component 2 (11.65%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                ),
                zaxis = list(
                    title = list(text = "<br><br>Component 3 (7.69%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                )
            ),
            legend = list(
                title = list(text = "Sample Type", 
                             font = list(size = 22, 
                                         face = "bold")),
                font = list(size = 18, 
                            face = "bold")
            )
        )
    
    ord_3d_plot
    
    # Save the plot as HTML
    saveWidget(ord_3d_plot, 
               file = "./Plots/KO_PCoA_By_SampleType_3D_Ordination_Plot.html", 
               selfcontained = TRUE)
    
    
    # save 3D plot
    # Load required library
    library(htmlwidgets)
    library(plotly)
    
    # Create the sample type plot
    ord_3d_plot <- plot_ly(
        data = Merged_df_weighted,
        x = ~ `Component 1 (19.9%)`,
        y = ~ `Component 2 (11.65%)`,
        z = ~ `Component 3 (7.69%)`,
        color = ~ SampleTypev2,
        colors = SampleType_custom_colors,
        text = ~ paste("SampleID:", SampleID,
                       "<br>RT Category:", RT_category,
                       "<br>Sample Type:", SampleTypev2),
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 4, alpha = 0.5)
    ) %>%
        layout(
            title = list(
                text = "Bray-Curtis PCoA based on KO abundances",
                font = list(size = 22, 
                            family = "Arial", 
                            color = "black")
            ),
            margin = list(t = 80), 
            scene = list(
                xaxis = list(
                    title = list(text = "<br><br>Component 1 (19.9%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                ),
                yaxis = list(
                    title = list(text = "<br><br>Component 2 (11.65%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                ),
                zaxis = list(
                    title = list(text = "<br><br>Component 3 (7.69%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                )
            ),
            legend = list(
                title = list(text = "Sample Type", 
                             font = list(size = 22, 
                                         face = "bold")),
                font = list(size = 18, 
                            face = "bold")
            )
        )
    
    ord_3d_plot
    
    # Save the plot as HTML
    saveWidget(ord_3d_plot, 
               file = "./Plots/KO_PCoA_By_SampleType_3D_Ordination_Plot.html", 
               selfcontained = TRUE)
    
    # Create the sample type plot
    ord_3d_plot <- plot_ly(
        data = Merged_df_weighted,
        x = ~ `Component 1 (19.9%)`,
        y = ~ `Component 2 (11.65%)`,
        z = ~ `Component 3 (7.69%)`,
        color = ~ RT_category,
        text = ~ paste("SampleID:", SampleID,
                       "<br>RT Category:", RT_category,
                       "<br>Sample Type:", SampleTypev2),
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 4, alpha = 0.5)
    ) %>%
        layout(
            title = list(
                text = "Bray-Curtis PCoA based on KO abundances",
                font = list(size = 22, 
                            family = "Arial", 
                            color = "black")
            ),
            margin = list(t = 80), 
            scene = list(
                xaxis = list(
                    title = list(text = "<br><br>Component 1 (19.9%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                ),
                yaxis = list(
                    title = list(text = "<br><br>Component 2 (11.65%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                ),
                zaxis = list(
                    title = list(text = "<br><br>Component 3 (7.69%)",
                                 standoff = 20,
                                 font = list(size = 18, 
                                             family = "Arial", 
                                             color = "black")),
                    tickfont = list(size = 16)
                )
            ),
            legend = list(
                title = list(text = "Sample Type", 
                             font = list(size = 22, 
                                         face = "bold")),
                font = list(size = 18, 
                            face = "bold")
            )
        )
    
    ord_3d_plot
    
    # Save the plot as HTML
    saveWidget(ord_3d_plot, 
               file = "./Plots/KO_PCoA_By_RT-cat_3D_Ordination_Plot.html", 
               selfcontained = TRUE)
    