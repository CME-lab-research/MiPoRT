# Script to join summary of PA (Presence-Absence) matrix with summary of abundance matrix for major groups.

library(tidyverse)

# Load the helper function
source("../Scripts/R_Utility_scripts/format_species_label.R")  # adjust path if needed

# Read presence-absence (P/A) summary file
PA_df_summary <- read.table("../05_a_Prevalence_analysis/Prevalence_Summary_filtered_Taxa_All_major_groups_HealthWise.txt", header = T, sep = "\t")
    # format species col
    PA_df_summary$Species <- str_replace_all(PA_df_summary$Species,
                                             "s__",
                                             "")

# Read species abundance summary file
Abundance_df_summary <- read.table("../05_a_Abundance_plots/Abundance_Taxa_filtered_min10_samples_All_MajorGroups_HealthWise.txt", header = T, sep = "\t")

# Abundance and prevlance summary was filtered for species present in min 10 samples, therefore has 780 species. Let's intersect both the tables for the same species.

# Check for common columns
intersect(colnames(Abundance_df_summary), colnames(PA_df_summary))

# Use left_join() to return all rows from x with a match in y.
Merged_summary_df <- left_join(PA_df_summary, Abundance_df_summary, join_by(Species == Species))

# sanity check
intersect(Abundance_df_summary$Species, PA_df_summary$Species)

# write this main merged table
write.table(Merged_summary_df, "./Summary_PA_and_Abundance_Species_min10.txt", 
            sep = "\t", 
            row.names = F,
            quote = F)

View(Merged_summary_df %>% select(c(Total_Prevalence.x, Total_Prevalence.y)))
# First, let's plot for Healthy samples
    Healthy_df <- Merged_summary_df %>%
        select(contains(c("Species", "Healthy_Abundance", "Healthy")))
    
    colnames(Healthy_df) <- c("Species", "Species_names", "Abundance_Upper.RT", "Abundance_Intermediate.RT", "Abundance_Lower.RT", "Abundance_Sputum", "Abundance_BAL", "Abundance_Other", "Abundance_Nasal.swab", "Abundance_Saliva", "Abundance_Buccal.mucosa", "Abundance_Oral.swab", "Abundance_Tongue.dorsum", "Prevalence_Upper.RT", "Prevalence_Intermediate.RT", "Prevalence_Lower.RT", "Prevalence_Nasal.swab", "Prevalence_Buccal.mucosa", "Prevalence_Oral.swab", "Prevalence_Saliva", "Prevalence_Tongue.dorsum", "Prevalence_Sputum", "Prevalence_BAL", "Prevalence_Other", "Prevalence_H"
)
    # format species names
    Healthy_df$Species_names <- format_species_label(Healthy_df$Species)
    
    # Melt all RT columns
    Healthy_df_long <- Healthy_df %>%
        pivot_longer(contains(c("Abundance", "Prevalence")) & ends_with("RT"), names_to = c("Metric", "RT_category"),  names_sep = "_", values_to = "Value") %>%
        pivot_wider(
            names_from = "Metric",
            values_from = "Value"
        ) %>%
        mutate(RT_category = factor(RT_category, 
                                       levels = c("Upper.RT", "Intermediate.RT", "Lower.RT" ),
                                       labels = c("Upper RT", "Intermediate RT",  "Lower RT" ))) %>%
        select(all_of(c("Species", "Species_names", "RT_category", "Abundance", "Prevalence")))

    # Get species order for plotting
    Species_Order <- Healthy_df_long %>%
        group_by(Species_names) %>%
        summarise(Avg = mean(Abundance),
                  Var_Abundance = var(Abundance)) %>%
        filter(Avg > 0.01) %>%
        arrange(desc(Var_Abundance)) %>%
        ungroup() %>%
        pull(Species_names) %>% head(30)
    
    Species_Order <- droplevels(Species_Order)
    
    Healthy_df_long_sub <- Healthy_df_long %>%
        filter(Species_names %in% Species_Order) %>%
        mutate(Species_names = factor(Species_names,
               levels = Species_Order,
               ordered = T))
"    
    Healthy_df_long$Species_names <- factor(Healthy_df_long$Species_names,
                                            levels = Species_Order)"
    
    # Horizontal bar plot for total prevalence for all Taxa
    ggObj <- ggplot(Healthy_df_long_sub, 
        aes(x = RT_category, 
            y = reorder(Species_names, Abundance), 
            fill = Abundance))
    
    savePlot <- ggObj + 
        geom_tile() + 
        labs(title="Mean abundance across healthy RT", x= "RT category", y= "Species", caption = "#Species: Top 60 species. Filter (Present >10 samples & \nMean abundance >0.1") + 
        scale_y_discrete(labels = function(x) parse(text = x)) +
        theme_bw() +
        theme(axis.text=element_text(size=22), 
              axis.text.x=element_text(size=22, angle = 45, vjust = 1, hjust = 1), 
              axis.title=element_text(size=22, face="bold"),
              title = element_text(size=24,face="bold"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 20, face = "bold"),
              #plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
              strip.text = element_text(size = 18, face = "bold")) +  
        scale_fill_gradientn(
            colours = c("#b3b6b7", "#1f78b4", "#eff088", "#feb24c", "#f03b20"),
            values = scales::rescale(c(0, 3, 6, 10, 15)),
            limits = c(0, 15),
            name = "Rel. Abundance",
            na.value = "grey"
        ) + coord_fixed()
    
    
    # save plot
    ggsave("./Plots/Total_abundance_Healthy_Species_min10_samples.png", 
           plot = savePlot,
           width = 50,              # Reduce dimensions
           height = 225,
           units = "cm",
           dpi = 300,               # High resolution
           limitsize = FALSE
    )
    

    # save plot with less Taxa
    ggsave("./Plots/Total_abundance_Healthy_Species_min10_samples_Top50.png", 
           plot = savePlot,
           width = 50,              # Reduce dimensions
           height = 125,
           units = "cm",
           dpi = 300,               # High resolution
           limitsize = FALSE
    )

    # Fig. 2. Bubble plot for prevalence and abundance

    # Make a newlist
    New_order <- reorder(Healthy_df_long_sub$Species_names, 
        Healthy_df_long_sub$Prevalence)
    
    # Use it to sort species on y
    ggObj <- ggplot(Healthy_df_long_sub, 
        aes(x = Prevalence, 
            y = New_order, 
            size = Abundance))
    
    #ggObj <- ggplot(Healthy_df_long_sub, aes(x = Prevalence, y =  Species_names, size = Abundance))
    
    
    savePlot <- ggObj + 
    geom_point(alpha = 0.7, color = "#009E73") + 
    facet_wrap(.~ RT_category) +
    scale_size_continuous(range = c(3, 10)) +  # Adjust bubble size range
    labs(
            title = "Species prevalence and abundance in healthy RT",
            x = "Prevalence (%)",
            y = "Species",
            size = "Mean Rel. abundance"
        ) +
    scale_y_discrete(labels = function(x) parse(text = x)) +
    theme_bw() +
    theme(axis.text=element_text(size=22), 
          axis.text.x=element_text(size=22, angle = 45, vjust = 1, hjust = 1), 
          axis.title=element_text(size=22, face="bold"),
          title = element_text(size=24,face="bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          #plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold")) 

    # save plot with less Taxa
    ggsave("./Plots/Bubble_plot_Abundance_Prevalence_Healthy_Top30_sp.png", 
           plot = savePlot,
           width = 50,              # Reduce dimensions
           height = 30,
           units = "cm",
           dpi = 300,               # High resolution
           limitsize = FALSE
    )

    
# Repeat for diseased samples
    Diseased_df <- Merged_summary_df %>%
        select(contains(c("Species", "Diseased_Abundance", "Diseased")))
    
    # sanity check
    colnames(Diseased_df)
    
    NewColNames <- c("Species", "Species_names", "Abundance_Upper.RT", "Abundance_Intermediate.RT", "Abundance_Lower.RT", "Abundance_Sputum", "Abundance_BAL", "Abundance_Supraglottal", "Abundance_Nasal.swab", "Abundance_Saliva", "Abundance_Buccal.mucosa", "Abundance_Oral.swab", "Prevalence_Upper.RT", "Prevalence_Intermediate.RT", "Prevalence_Lower.RT", "Prevalence_Nasal.swab", "Prevalence_Buccal.mucosa", "Prevalence_Oral.swab", "Prevalence_Saliva", "Prevalence_Supraglottal", "Prevalence_Sputum", "Prevalence_BAL", "Prevalence_Other", "Prevalence_D")

    colnames(Diseased_df) <- NewColNames
        
    # format species names
    Diseased_df$Species_names <- format_species_label(Diseased_df$Species)
    
    # Melt all RT columns
    Diseased_df_long <- Diseased_df %>%
        pivot_longer(contains(c("Abundance", "Prevalence")) & ends_with("RT"), names_to = c("Metric", "RT_category"),  names_sep = "_", values_to = "Value") %>%
        pivot_wider(
            names_from = "Metric",
            values_from = "Value"
        ) %>%
        mutate(RT_category = factor(RT_category, 
                                    levels = c("Upper.RT", "Intermediate.RT", "Lower.RT" ),
                                    labels = c("Upper RT", "Intermediate RT",  "Lower RT" ))) %>%
        select(all_of(c("Species", "Species_names", "RT_category", "Abundance", "Prevalence")))
    
    # Get species order for plotting
    Species_Order <- Diseased_df_long %>%
        group_by(Species_names) %>%
        summarise(Avg = mean(Abundance),
                  Var_Abundance = var(Prevalence),
                  Prevalence = Prevalence) %>%
        filter(Avg > 0.01) %>%
        arrange(desc(Prevalence)) %>%
        ungroup() %>%
        pull(Species_names) %>% head(30)
    
    Diseased_df_long_sub <- Diseased_df_long %>%
        filter(Species_names %in% Species_Order) %>%
        mutate(Species_names = factor(Species_names,
                                      levels = Species_Order,
                                      ordered = T))
    # Fig. 2. Bubble plot for prevalence and abundance
    
    #ggObj <- ggplot(Diseased_df_long_sub, aes(x = Prevalence, y = reorder(Species_names, Prevalence), size = Abundance))
    
    # Add custom sp. order list
    
    ggObj <- ggplot(Diseased_df_long_sub, 
        aes(x = Prevalence, 
            y = New_order, 
            size = Abundance))
    
    savePlot <- ggObj + 
        geom_point(alpha = 0.7, color = "#CC79A7") +
        facet_wrap(.~ RT_category) +
        scale_size_continuous(range = c(3, 10)) +  # Adjust bubble size range
        labs(
            title = "Species prevalence and abundance in diseased RT",
            x = "Prevalence (%)",
            y = "Species",
            size = "Mean Rel. abundance"
        ) +
        scale_y_discrete(labels = function(x) parse(text = x)) +
        theme_bw() +
        theme(axis.text=element_text(size=22), 
              axis.text.x=element_text(size=22, angle = 45, vjust = 1, hjust = 1), 
              axis.title=element_text(size=22, face="bold"),
              title = element_text(size=24,face="bold"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 20, face = "bold"),
              #plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
              strip.text = element_text(size = 18, face = "bold")) 
    
    # save plot with less Taxa
    ggsave("./Plots/Bubble_plot_Abundance_Prevalence_Diseased_Top30_sp.png", 
           plot = savePlot,
           width = 50,              # Reduce dimensions
           height = 30,
           units = "cm",
           dpi = 300,               # High resolution
           limitsize = FALSE
    )

    
    # can we merge those lists? plots?
colnames(Healthy_df_long_sub)  
[1] "Species"       "Species_names" "RT_category"   "Abundance"     "Prevalence"   
[6] "Health_Status"

colnames(Diseased_df_long_sub)
[1] "Species"       "Species_names" "RT_category"   "Abundance"     "Prevalence"   
[6] "Health_Status"

Healthy_df_long_sub$Health_Status <- rep("Healthy", 90)
Diseased_df_long_sub$Health_Status <- rep("Diseased", 90)

Combined_df <- bind_rows(Healthy_df_long_sub, Diseased_df_long_sub)

# left join these 2
Combined_df$Health_Status <- factor(Combined_df$Health_Status,
                                    levels = c("Healthy", "Diseased"))

ggObj <- ggplot(Combined_df, aes(x = Prevalence, y = reorder(Species_names,  Prevalence), size = Abundance, colour = Health_Status))

Manual_colors <- c("Healthy" = "#009E73", 
                   "Diseased"= "#CC79A7")

savePlot <- ggObj + 
    geom_point(alpha = 0.7) +
    facet_wrap(.~ RT_category) +
    scale_size_continuous(range = c(3, 10)) +  # Adjust bubble size range
    labs(
        title = "Species prevalence and abundance across RT",
        x = "Prevalence (%)",
        y = "Species",
        size = "Mean Rel. abundance"
    ) +
    scale_colour_manual(values = Manual_colors) +
    scale_y_discrete(labels = function(x) parse(text = x)) +
    theme_bw() +
    theme(axis.text=element_text(size=22), 
          axis.text.x=element_text(size=22, angle = 45, vjust = 1, hjust = 1), 
          axis.title=element_text(size=22, face="bold"),
          title = element_text(size=24,face="bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          #plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold")) 

# save plot with less Taxa
ggsave("./Plots/Bubble_plot_Abundance_Prevalence_Top30_sp.png", 
       plot = savePlot,
       width = 50,              # Reduce dimensions
       height = 30,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE
)
    