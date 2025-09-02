# Script for upset plot from Prevalence feature table

# Load libraries
library(ComplexUpset)
library(tidyverse)

# Read presence absences table & metadata table
Metaphlan_PA_df <- read.table("../MetaPhlan4_results/MiPORT_filtered_min10_samples_PA_data.txt", header = T, sep = "\t")

# Move rownames as column
Metaphlan_PA_df <- Metaphlan_PA_df %>%
    rownames_to_column(var = 'Species')

Species_of_Interest <- row.names(Metaphlan_PA_df)

rm(Metaphlan_PA_df)

# Read raw PA table to add Anterior Nare samples
    #Raw_PA_df <- read.table("../MetaPhlan4_results/MiPORT_raw_PA_data.txt", sep = "\t", header = T)

# Read Metadata
Metadata_Df <- 
    read.table("../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", header = T, sep = "\t")

# Find samples present in both files to retain them
Retain_Samples <- 
    intersect(colnames(Metaphlan_PA_df), 
              Metadata_Df$SampleID)

# Filter Raw PA table for samples of interest
Raw_PA_df <- Metaphlan_PA_df %>%
    select(all_of(Retain_Samples))

Metadata_Df <- Metadata_Df %>%
    filter(SampleID %in% Retain_Samples)

table(Metadata_Df$SampleTypev2, useNA = 'ifany')

    # Use this table for PA with Anterior nares samples (Species are filtered for min 10 sample presence)
    write.table(Raw_PA_df, "../MetaPhlan4_results/MiPORT_filtered_withAN_PA_data.txt", row.names = T, quote = F, sep = "\t")

# Add factors to metadata
    # Factors for SampleTypev2
    SamplingSite_Factor <- c("Anterior_nares", "Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Sputum", "BAL", "Other")
    SamplingSite_labels <- c("Anterior nares", "Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL", "Other")
    
    # Add factors for SampleTypev2
    Metadata_Df$SampleTypev2 <- 
        factor(Metadata_Df$SampleTypev2, 
               levels = SamplingSite_Factor, 
               labels = SamplingSite_labels)  
    
    levels(Metadata_Df$SampleTypev2)

    # Add factors for RT category
    table(Metadata_Df$RT_category, useNA = 'ifany')
    RT_Cat_Factor <- c("URT", "IRT", "LRT")
    RT_Cat_Labels <- c("Upper RT", "Intermediate RT", "Lower RT")
    Metadata_Df$RT_category <- factor(Metadata_Df$RT_category , 
                                      levels = RT_Cat_Factor,
                                      labels = RT_Cat_Labels)
    
    # sanity check
    table(Metadata_Df$RT_category, useNA = 'ifany')
    
    # Add factors for health status
    table(Metadata_Df$Healthy, useNA = 'ifany')

    Health_stat_lvls <- c("TRUE", "FALSE", "Unknown")
    Health_stat_labs <- c("Healthy", "Diseased", "Unknown")
    
    Metadata_Df$Healthy <- factor(Metadata_Df$Healthy, 
                                  levels = Health_stat_lvls,
                                  labels = Health_stat_labs)
    # sanity check
    table(Metadata_Df$Healthy, useNA = 'ifany')
    
    Metadata_Df_sub <- Metadata_Df %>%
        select(all_of(c("SampleID", "Healthy", "RT_category", "SampleTypev2")))
    
# Step 1. Convert wide P/A table into long format
    Species_Presence_long <- Metaphlan_PA_df %>%
        pivot_longer(
            cols = -Species, # All columns except Species
            names_to = "SampleID",
            values_to = "Present"
        ) %>%
        filter(Present == 1)  # Keep only present species
    
    # Add metadata information to this table
    Merged_df <- left_join(Species_Presence_long, Metadata_Df_sub,
                           by = join_by("SampleID" == "SampleID"))

# Step 2. For each species, collapse into list of RT-categories it's found in

    # Mark a species as present in a RT category only if it’s seen in, ≥10% of samples in that RT-category
    Species_PA_RTCAT <- Merged_df %>%
        group_by(Species, RT_category) %>%
        summarise(prev_RT = mean(Present), .groups = "drop") %>%
        mutate(Present_RT = as.integer(prev_RT >= 0.10)) %>%
        select(-prev_RT) %>%
        pivot_wider(names_from = RT_category, 
                    values_from = Present_RT, 
                    values_fill = 0) 

    # Make the same for sampleType variable and merge 2 tables
    Species_PA_SampleType <- Merged_df %>%
        group_by(Species, SampleTypev2) %>% 
        summarise(prev_ST2 = mean(Present), .groups = "drop") %>%
        mutate(Present_ST2 = as.integer(prev_ST2 >= 0.10)) %>%
        select(-prev_ST2) %>%
        pivot_wider(names_from = SampleTypev2, 
                    values_from = Present_ST2, 
                    values_fill = 0) 
    
# Merge and upset plot RT category summary
Plot_df <- left_join(Species_PA_RTCAT, 
                     Species_PA_SampleType, by = "Species")

    # sanity check
    sapply(Plot_df[RT_Cat_Labels], function(x) all(x %in% c(0,1,TRUE,FALSE)))
    # All mus tbe binary values so TRUE

## SAVE THIS plot df for plotting Upsetfaster next time
write.table(Plot_df,
            "Upset_plot_matrix_Species_Present_whenPrevalence_gt10_percent.txt",
            row.names = F,
            sep = '\t')


# Add variable names to plot in a list
Cols_to_Plot <- c(RT_Cat_Labels, SamplingSite_labels)
RT_Cat_order <- list(
    c("Upper RT", "Intermediate RT", "Lower RT"),
    c("Upper RT", "Lower RT"),
    c("Intermediate RT", "Lower RT"),
    c("Upper RT", "Intermediate RT"),
    "Upper RT",
    "Intermediate RT",
    "Lower RT"
)

# Plot 1 RT_category
Plot_Obj <- upset(
    Plot_df,
    intersect = RT_Cat_Labels,            # the set columns
    intersections = RT_Cat_order,         # ordered list of combos
    name = "RT-category intersection",
    width_ratio = 0.3,
    height_ratio= 0.4,
    sort_intersections = FALSE,
    sort_sets = 'descending',
    warn_when_dropping_groups = TRUE,
    min_size = 0,                         # allow empty if you requested them
    themes = upset_default_themes(
        text = ggplot2::element_text(size = 22, face = "bold")
    )
) 

Plot_Obj

# Add custom theme
savePlot <- Plot_Obj +
    theme_bw() +
    theme(
        axis.text.y    = element_text(size = 22),
        axis.text.x    = element_blank(),
        axis.title.y   = element_blank(),
        title          = element_text(size = 24, face = "bold"),
        legend.text    = element_text(size = 20),
        legend.title   = element_text(size = 20, face = "bold"),
        plot.margin    = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text     = element_text(size = 18, face = "bold")
    )


savePlot

# save png plot
    ggsave("./Plots/Upset_RT_categories_with_species_min10_samples.png", 
           plot = savePlot,
           width = 30,              # Reduce dimensions
           height = 20,
           units = "cm",
           dpi = 600,               # High resolution
           limitsize = FALSE
    )

# save svg
ggsave("./Plots/Upset_RT_categories_with_species_min10_samples.svg", 
       plot = savePlot,
       width = 30,              # Reduce dimensions
       height = 20,
       units = "cm",
       limitsize = FALSE
)

# Plot 2 all Sample type species intersections
# Without Other & AN
SamplingSite_labels_sub <- c("Anterior nares", "Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL")

Plot_df_sub <- Plot_df[,SamplingSite_labels_sub]
Species_to_rm <- (rowSums(Plot_df_sub) != 0)

table(Species_to_rm)

Plot_df_sub <- Plot_df_sub[Species_to_rm,]

Plot_Obj_2 <- upset(
    Plot_df_sub,
    intersect = SamplingSite_labels_sub,            # the set columns
    name = "Sample-type intersection",
    width_ratio = 0.2,
    height_ratio= 0.5,
    sort_intersections = 'descending',
    sort_sets = FALSE,
    warn_when_dropping_groups = TRUE,
    min_size = 0,                         # allow empty if you requested them
    themes = upset_default_themes(
        text = ggplot2::element_text(size = 22, face = "bold")
    )
) 

    Plot_Obj_2
    # Add themes
    savePlot_2 <- Plot_Obj_2 +
        theme_bw() +
        theme(
            axis.text.y    = element_text(size = 22),
            axis.text.x    = element_blank(),
            axis.title.y   = element_blank(),
            title          = element_text(size = 24, face = "bold"),
            legend.text    = element_text(size = 20),
            legend.title   = element_text(size = 20, face = "bold"),
            plot.margin    = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
            strip.text     = element_text(size = 18, face = "bold")
        )

    # Check plot
    savePlot_2

# save plot
ggsave("./Plots/Global_Upset_Sampletypes_with_531species_min10_samples.png", 
       plot = savePlot_2,
       width = 45,              # Reduce dimensions
       height = 25,
       units = "cm",
       dpi = 600,               # High resolution
       limitsize = FALSE
)

# save svg
ggsave("./Plots/Global_Upset_Sampletypes_with_531species_min10_samples.svg", 
       plot = savePlot_2,
       width = 45,              # Reduce dimensions
       height = 25,
       units = "cm",
       limitsize = FALSE
)

# Plot sampletype plot with filtered annotations
Plot_Obj_2 <- upset(
    Plot_df,
    intersect = SamplingSite_labels,            # the set columns
    name = "Sample-type intersection",
    width_ratio = 0.2,
    height_ratio= 0.5,
    sort_intersections = 'descending',
    sort_sets = FALSE,
    warn_when_dropping_groups = TRUE,
    min_size = 0,                         # allow empty if you requested them
    themes = upset_default_themes(
        text = ggplot2::element_text(size = 22, face = "bold")
    )
) 

Plot_Obj_2
# Add themes
savePlot_2 <- Plot_Obj_2 +
    theme_bw() +
    theme(
        axis.text.y    = element_text(size = 22),
        axis.text.x    = element_blank(),
        axis.title.y   = element_blank(),
        title          = element_text(size = 24, face = "bold"),
        legend.text    = element_text(size = 20),
        legend.title   = element_text(size = 20, face = "bold"),
        plot.margin    = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text     = element_text(size = 18, face = "bold")
    )

# Check plot
savePlot_2

# save plot
ggsave("./Plots/Global_Upset_Sampletypes_with_species_min10_samples.png", 
       plot = savePlot_2,
       width = 45,              # Reduce dimensions
       height = 25,
       units = "cm",
       dpi = 600,               # High resolution
       limitsize = FALSE
)

# save svg
ggsave("./Plots/Global_Upset_Sampletypes_with_species_min10_samples.svg", 
       plot = savePlot_2,
       width = 45,              # Reduce dimensions
       height = 25,
       units = "cm",
       limitsize = FALSE
)


Plot_df %>%
    filter(rowSums(across(SamplingSite_labels)) == (length(SamplingSite_labels) -2))

# all combinations of all subset lengths
all_combinations <- function(x, min_k = 1, max_k = length(x)) {
    x <- unique(as.character(x))  # ensure uniqueness
    out <- list()
    for (k in seq(min_k, max_k)) {
        cmbs <- combn(x, k, simplify = FALSE)
        out <- c(out, cmbs)
    }
    out
}

Sampletype_order <- all_combinations(SamplingSite_labels)

# define custom intersection combinations
Sampletype_order <- list(
    # all
    c("Anterior nares", "Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL"),
    
    # w/o AN
    c("Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL"),
    
    # without NS
    c("Anterior nares", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL"),
    c("Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL"),
    c("Anterior nares", "Supraglottal"),
    
    # URT
    c("Buccal mucosa", "Oral swab", "Saliva"),
    
    # URT-IRT
    c("Buccal mucosa", "Saliva", "Tongue dorsum"),
    c("Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum"),
    c("Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Sputum"),
    c("Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal"),
    
    # URT-LRT
    c("Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Sputum", "BAL"),
    
    # each sample types
    "Anterior nares", "Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL"
)


# Plot sampletype plot with filtered annotations
Plot_Obj_2 <- upset(
    Plot_df_sub,
    intersect = SamplingSite_labels_sub, # the set columns
    intersections = Sampletype_order, 
    name = "Sample-type intersection",
    width_ratio = 0.3,
    height_ratio= 0.5,
    sort_intersections = FALSE,
    sort_sets = FALSE,
    warn_when_dropping_groups = FALSE,
    min_size = 0,                         # allow empty if you requested them
    themes = upset_default_themes(
        text = ggplot2::element_text(size = 22, face = "bold")
    )
) 

Plot_Obj_2

# Add themes
savePlot_2 <- Plot_Obj_2 +
    theme_bw() +
    theme(
        axis.text.y    = element_text(size = 22),
        axis.text.x    = element_blank(),
        axis.title.y   = element_blank(),
        title          = element_text(size = 24, face = "bold"),
        legend.text    = element_text(size = 20),
        legend.title   = element_text(size = 20, face = "bold"),
        plot.margin    = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text     = element_text(size = 18, face = "bold")
    )

# Check plot
savePlot_2

# save plot
ggsave("./Plots/Upset_Sampletypes_with_species_min10_samples.png", 
       plot = savePlot_2,
       width = 40,              # Reduce dimensions
       height = 25,
       units = "cm",
       dpi = 600,               # High resolution
       limitsize = FALSE
)

# save svg
ggsave("./Plots/Upset_Sampletypes_with_species_min10_samples.svg", 
       plot = savePlot_2,
       width = 35,              # Reduce dimensions
       height = 25,
       units = "cm",
       limitsize = FALSE
)


# Without Other & AN
SamplingSite_labels_sub <- c("Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL")
Plot_df_sub <- Plot_df[,SamplingSite_labels_sub]
    
Species_Rs <- rowSums(Plot_df_sub) !=0

# Plot sampletype plot with filtered annotations
Plot_Obj_2 <- upset(
    Plot_df_sub[Species_Rs,],
    intersect = SamplingSite_labels_sub, # the set columns
    #intersections = Sampletype_order, 
    name = "Sample-type intersection",
    width_ratio = 0.3,
    height_ratio= 0.5,
    sort_intersections = 'descending',
    sort_sets = FALSE,
    warn_when_dropping_groups = FALSE,
    min_size = 0,                         # allow empty if you requested them
    themes = upset_default_themes(
        text = ggplot2::element_text(size = 22, face = "bold")
    )
) 

Plot_Obj_2

# Add themes
savePlot_2 <- Plot_Obj_2 +
    theme_bw() +
    theme(
        axis.text.y    = element_text(size = 22),
        axis.text.x    = element_blank(),
        axis.title.y   = element_blank(),
        title          = element_text(size = 24, face = "bold"),
        legend.text    = element_text(size = 20),
        legend.title   = element_text(size = 20, face = "bold"),
        plot.margin    = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text     = element_text(size = 18, face = "bold")
    )

# Check plot
savePlot_2

# save plot
ggsave("./Plots/v2_Upset_Sampletypes_with_species_min10_samples.png", 
       plot = savePlot_2,
       width = 40,              # Reduce dimensions
       height = 25,
       units = "cm",
       dpi = 600,               # High resolution
       limitsize = FALSE
)

# save svg
ggsave("./Plots/v2_Upset_Sampletypes_with_species_min10_samples.svg", 
       plot = savePlot_2,
       width = 35,              # Reduce dimensions
       height = 25,
       units = "cm",
       limitsize = FALSE
)



# Map sample type → RT category
sampletype_to_rt <- c(
    "Anterior nares" = "Upper RT",
    "Nasal swab"     = "Upper RT",
    "Buccal mucosa"  = "Upper RT",
    "Oral swab"      = "Upper RT",
    "Saliva"         = "Upper RT",
    "Tongue dorsum"  = "Intermediate RT",
    "Supraglottal"   = "Intermediate RT",
    "Sputum"         = "Intermediate RT",
    "BAL"            = "Lower RT"
)








