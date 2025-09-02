# Script to do merge alpha diversity and metadata columns to do some statistics
library(tidyverse)
library(introdataviz)
library(ggpubr)


# Part 1: Merge metadata and alpha diversity metric values

  # read Metadata table
  Metadata_df <- read.table("../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", header = T, sep = "\t")

    # remove blank rows
    Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]
    Alpha_divs <- c("Gini-Dominance", "Richness", "Shannon", "Gini-Simpson")

    # Rm existing alpha div values
    Metadata_df <- Metadata_df %>% 
      #select(!Alpha_divs) %>%
      relocate(SampleTypev2, .after = 'SampleType')

      colnames(Metadata_df)

# List all diversity file names
  Alpha_div_files <- c("../MetaPhlan4_results/Diversity_merged_Batch_corrected_data/MiPORT_All_RT_featureTable_sp_filt_1e-4_Relab_per_Min10_batchCorrected_samples_gini.tsv", 
                       "../MetaPhlan4_results/Diversity_merged_Batch_corrected_data/MiPORT_All_RT_featureTable_sp_filt_1e-4_Relab_per_Min10_batchCorrected_samples_richness.tsv",
                       "../MetaPhlan4_results/Diversity_merged_Batch_corrected_data/MiPORT_All_RT_featureTable_sp_filt_1e-4_Relab_per_Min10_batchCorrected_samples_shannon.tsv ",
                       "../MetaPhlan4_results/Diversity_merged_Batch_corrected_data/MiPORT_All_RT_featureTable_sp_filt_1e-4_Relab_per_Min10_batchCorrected_samples_simpson.tsv")

  Merged_df <- Metadata_df

  # read each div file and add values to merged df
  for(eachDiv in Alpha_div_files){
    print(paste0("Processing file: ", eachDiv))
    # Reading diversity file
    Diversity_df <- read.table(eachDiv, sep = '\t', header = T, row.names = 1)
    
    alphaDiv <- colnames(Diversity_df)
    
    Diversity_df <- Diversity_df %>% rownames_to_column(var = "SampleID")
    colnames(Diversity_df) <- c("SampleID", alphaDiv)
    
    print(alphaDiv)
    
    # left join data matrix and sample IDs
    Merged_df <- Merged_df %>%
      left_join(Diversity_df , by = join_by("SampleID" == "SampleID"))%>% 
      relocate(alphaDiv, .after = last_col())
    
    print(colnames(Merged_df))
  }

  # Rename columns
  colnames(Merged_df) <-
    c(colnames(Merged_df)[1:28], 'Gini-Dominance', 'Richness',
      'Shannon', 'Gini-Simpson')
    # sanity check
    colnames(Merged_df)

  # Check if you have QC failed samples or samples with failed metaphlan profiling
  Failed_Samples <- Merged_df %>% 
    select(all_of(Alpha_divs)) %>% 
    apply(MARGIN = 1, FUN = function(eachRow){
    sum(is.na(eachRow))
  })

  # sanity check
  table(Failed_Samples)
  '   0 
    3135
  No samples have missing/NA values. Exclude 10 Unknown samples. 
  '
  # Update the same metadata table with alpha diversity values
  write.table(Merged_df, 
              "../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", 
              row.names = F, sep = '\t', quote = F)

# Remove AN samples
Merged_df_sub <- Merged_df %>% filter(SampleType != "Anterior_nares")

# SKIP: Add a new sampletype col which only includes top 8 sampletypes
  # Count frequency of each SampleType
  top_sampletypes <- Merged_df_sub %>%
    count(SampleType, sort = TRUE) %>%   # Count and sort
    slice_max(n, n = 8) %>%              # Get top 8 most frequent
    pull(SampleType)                     # Extract SampleType names
  
    # Check
    top_sampletypes

  # Re-assign remaining sampletypes as "Other"
  Merged_df <- Merged_df %>%
    mutate(SampleTypev2 = ifelse(as.character(SampleType) %in% top_sampletypes, as.character(SampleType), "Other"))
  # Check dims
  dim(Merged_df_sub)
  
  # Make a long df for plotting diversity distribution
  Merged_df_long <- Merged_df_sub %>% 
    pivot_longer(cols = Alpha_divs, 
                 names_to = "Alpha_div", 
                 values_to = "Div_Value")
  
  # Check dims
  dim(Merged_df_long)
  # should be 2927 * 4 > 11708
  
  table(Merged_df_sub$SampleTypev2, useNA = 'ifany')
    '
Anterior_nares            BAL  Buccal_mucosa     Nasal_Swab      Oral_swab 
           208            578           1027             27            246 
         Other         Saliva         Sputum   Supraglottal  Tongue_dorsum 
            64            315            193             59            418 
    '
 
# Add factors for SampleTypev2
  # Make an ordered list for factors
  SamplingSite_Factor <- c("Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Sputum", "BAL", "Other")

  # Format factor text labels
  SampleTypev2_labels <- c(
    #"Anterior nares\n(n=208)",
    "Nasal swab\n(n=27)",
    "Buccal mucosa\n(n=1027)",
    "Oral swab\n(n=246)",
    "Saliva\n(n=315)",
    "Tongue dorsum\n(n=418)",
    "Supraglottal\n(n=59)",
    "Sputum\n(n=193)",
    "BAL\n(n=578)",
    "Other\n(n=64)"
  )
  
  SampleTypev2_labels_2 <- c(
    #"Anterior nares\n(n=208)",
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

  # For plotting define custom colors for top 9 + "Other"
  SampleType_custom_colors <- c(
    #"Anterior nares\n(n=208)" = "#00BFC4",
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
  
  # Add factors
  Merged_df_long$SampleTypev2 <- factor(Merged_df_long$SampleTypev2, 
                                   levels = SamplingSite_Factor,
                                   labels = SampleTypev2_labels)  

  levels(Merged_df_long$SampleTypev2)

  Merged_df_long$SampleTypev2 <- factor(Merged_df_long$SampleTypev2, 
                                        levels = SampleTypev2_labels,
                                        labels = SampleTypev2_labels_2)  
  # sanity check
  table(Merged_df_long$SampleTypev2, useNA = 'ifany')/4

# Add factors for RT category
  table(Merged_df_long$RT_category, useNA = 'ifany')/4
  '
   IRT  LRT  URT 
   719  578 1630
  '
  # specify factor vars
  RT_Cat_Factor <- c("URT", "IRT", "LRT")
  RT_Cat_Labels <- c("Upper RT (N=1630)", 
                     "Intermediate RT (N=719)", 
                     "Lower RT (N=578)")
  
  Merged_df_long$RT_category <- factor(Merged_df_long$RT_category , 
                                 levels = RT_Cat_Factor,
                                 labels = RT_Cat_Labels)

  # sanity check
  table(Merged_df_long$RT_category, useNA = 'ifany')/4

# Add factors for health status
  table(Merged_df_long$Healthy, useNA = 'ifany')/4
  
  Health_stat_lvls <- c("TRUE", "FALSE", "Unknown")
  Health_stat_labs <- c("Healthy (N=1460)", 
                        "Diseased (N=1457)", 
                        "Unknown (N=10)")
  
  Merged_df_long$Healthy <- factor(Merged_df_long$Healthy, 
                               levels = Health_stat_lvls,
                               labels = Health_stat_labs)
  # sanity check
  table(Merged_df_long$Healthy, useNA = 'ifany')/4

  # set colors for plots
  Grp_Healthy_custom_color <- c(
    "Healthy (N=1460)" = "#59D7D9",  # Sky blue
    "Diseased (N=1457)" = "#E6550D"  # A rich burnt orange
  )
  
# filter out Unknown health status samples
  Merged_df_long <- Merged_df_long %>% filter(Healthy != "Unknown (N=10)")
  Merged_df_long$Healthy <- droplevels(Merged_df_long$Healthy)

# Plot 1. Alpha diversity density
ggplotObj <- ggplot(Merged_df_long, 
                    aes(Div_Value, 
                        fill = Healthy))

savePlot <- ggplotObj + 
  geom_density(kernel = "gaussian", alpha = 0.8) + 
  facet_wrap(Alpha_div ~ ., scales = "free", nrow = 2) +
  labs(title = "Alpha diversity distribution for health status", 
       fill = "Health status (N)",
       y = "Density",
       x = "Alpha diversity") + 
  scale_fill_manual(values = Grp_Healthy_custom_color) +
  theme_bw() +
  theme(axis.text=element_text(size=22), 
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=24,face="bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text = element_text(size = 18, face = "bold"),
        strip.background = element_blank())  

  savePlot
  
  # save png plot
    ggsave("./Batch_corrected_Plots/Alpha_div_dist_byHealth.png", 
           plot = savePlot,
           width = 32,              # Reduce dimensions
           height = 20,
           units = "cm",
           limitsize = FALSE,
           dpi = 300               # High resolution
    )

  # save svg plot
    ggsave("./Batch_corrected_Plots/Alpha_div_dist_byHealth.svg", 
           plot = savePlot,
           width = 32,              # Reduce dimensions
           height = 20,
           units = "cm",
           limitsize = FALSE
    )

# Plot 2 freq RT cat wrt Health
  # Custom labels for facets
  Alpha_div_custom_labels <- c(
    "Gini.dominance" = "Gini_dominance",
    "Richness" = "Richness",
    "Shannon" = "Shannon",
    "Gini.simpson" = "Gini_simpson"
  )

  # Plot obj
  ggplotObj <- ggplot(Merged_df_long, aes(Div_Value, colour = RT_category,
                                        linetype = Healthy))

  savePlot <- ggplotObj + 
    geom_freqpoly(linewidth = 1.5) + 
    facet_wrap(~ Alpha_div, 
               nrow = 2, ncol = 2, 
               scales = "free" 
               #labeller = labeller(Alpha_div = Alpha_div_custom_labels)
               ) +
    labs(
      title = "Alpha diversity distribution for all samples (N=2927)",
      colour = "RT category",
      x = "Diversity",
      y = "Sample Count"
    ) + 
    theme_bw() +
    theme(
      axis.text = element_text(size = 22), 
      axis.title = element_text(size = 22, face = "bold"),
      title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_blank()
    )

  savePlot

  # save png plot
  ggsave("./Batch_corrected_Plots/Alpha_div_dist_by_RT.png", 
         plot = savePlot,
         width = 35,              # Reduce dimensions
         height = 20,
         units = "cm",
         dpi = 300,               # High resolution
         limitsize = FALSE
  )

  # save svg plot
  ggsave("./Batch_corrected_Plots/Alpha_div_dist_by_RT.svg", 
         plot = savePlot,
         width = 35,              # Reduce dimensions
         height = 20,
         units = "cm",
         limitsize = FALSE
  )

# Plot 3 freq with sample types
  ggplotObj <- ggplot(Merged_df_long, aes(Div_Value, colour = SampleTypev2))

savePlot <- ggplotObj + geom_freqpoly(linewidth = 1.5) + 
  facet_wrap(Alpha_div ~ RT_category, 
             scales = "free", 
             nrow = 4) +
  scale_color_manual(values = SampleType_custom_colors) +
  labs(title = "Alpha diversity distribution for all samples (N=2927)",
       colour = "Sample type (N)",
       x = "Diversity",
       y = "Sample Count") + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 22), 
    axis.title = element_text(size = 22, face = "bold"),
    title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
    strip.text = element_text(size = 18, face = "bold"),
    strip.background.x = element_blank()
  )

savePlot

  # save png plot
  ggsave("./Batch_corrected_Plots/Alpha_div_dist_by_SampleTypev2.png",
         plot = savePlot,
         width = 50,          
         height = 30,
         units = "cm",
         dpi = 300,               # High resolution
         limitsize = FALSE)
  
  # save png plot
  ggsave("./Batch_corrected_Plots/Alpha_div_dist_by_SampleTypev2.svg",
         plot = savePlot,
         width = 50,          
         height = 30,
         units = "cm",
         limitsize = FALSE)

# Plot 4: Box plot diversity by RT category

  table(Merged_df_long$Alpha_div)
  '
  Gini-Dominance   Gini-Simpson       Richness        Shannon 
          2917           2917           2917           2917
  
  '
  
  table(Merged_df_long$RT_category)
  '
       Upper RT (N=1630) Intermediate RT (N=719)        Lower RT (N=578) 
                   6480                    2876                    2312 
  
  '

  # statistics

  # What to compare list?
  Comparisons <- list(
    c("Upper RT (N=1630)", "Intermediate RT (N=719)"),
    c("Upper RT (N=1630)", "Lower RT (N=578)"),
    c("Intermediate RT (N=719)", "Lower RT (N=578)")
  )
  
  ggplotObj <- ggplot(Merged_df_long, 
                      aes(x = RT_category,
                          y = Div_Value
                          ))
  
# save plot
  savePlot <- ggplotObj + geom_boxplot() +
    labs(title = "Alpha diversity across all respiratory tracts (N=2927)",
         colour = "QC status",
         x = NULL,
         y = "Diversity Value") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text = element_text(size = 22), 
          axis.title = element_text(size = 22, face = "bold"),
          title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold"),
          strip.background.x = element_blank()) + 
    facet_wrap(Alpha_div ~ ., scales = "free", nrow = 1) + 
    stat_compare_means(comparisons = Comparisons,
                       group.by = "Alpha_div",
                       label = "p.format",
                       method = "wilcox.test",
                       hide.ns = TRUE,
                       p.adjust.method = "BH",
                       show.legend = TRUE)

savePlot

  # save png plot
  ggsave("./Batch_corrected_Plots/Alpha_div_boxplots_wth_Stats_by_RT.png",
         plot = savePlot,
         width = 35,          
         height = 20,
         units = "cm",
         dpi = 300,               # High resolution
         limitsize = FALSE)

  # save svg plot
  ggsave("./Batch_corrected_Plots/Alpha_div_boxplots_wth_Stats_by_RT.svg",
         plot = savePlot,
         width = 35,          
         height = 20,
         units = "cm",
         limitsize = FALSE)
  
# Plot 5: Boxplot diversity by sampletypev2
table(Merged_df_sub$SampleTypev2)

# plot by RT cat
ggplotObj <- ggplot(Merged_df_long %>% 
                      filter(SampleTypev2 != 'Other (n=64)' ), 
                    aes(x = SampleTypev2, 
                        y = Div_Value, 
                        fill = SampleTypev2
                        ))

savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity for major sample types",
       fill = "Sample Type",
       x = NULL,
       y = "Diversity Value") + 
  scale_fill_manual(values = SampleType_custom_colors) +
  theme_bw() +
  facet_wrap(Alpha_div ~ ., scales = "free", nrow = 2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22, face = "bold"),
        title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        plot.margin = margin(t = 20, r = 20, b = 80, l = 20, unit = "pt"),
        strip.text = element_text(size = 18, face = "bold"),
        strip.background.x = element_blank()) 

savePlot

# save plot
ggsave("./Plots/Alpha_div_boxplots_by_sampletype_updated.png",
       width = 30,          
       height = 22,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE)

ggsave("./Plots/Alpha_div_boxplots_by_sampletype_updated.svg",
       width = 30,          
       height = 22,
       units = "cm")

# strip not working.
# Plot with strip
# Step 1: Create label strip data
# Unique mapping of sample type to RT category
strip_labels_df <- Merged_df_long %>%
  select(SampleTypev2, RT_category) %>%
  distinct() %>%
  mutate(x_pos = as.integer(SampleTypev2)) %>%
  group_by(RT_category) %>%
  summarise(x_start = min(x_pos),
            x_end = max(x_pos),
            label_pos = (x_start + x_end) / 2,
            .groups = "drop")

# Add a y position just below the plot (you can tweak)
strip_labels_df$x_start <- c(1, 7, 8.7)
strip_labels_df$x_end <- c(6.5, 8.5, 9.5)
strip_labels_df$label_pos <- c(2, 5, 8)


savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity for all samples",
       fill = "Sample Type",
       x = NULL,
       y = "Diversity Value") + 
  scale_fill_manual(values = SampleType_custom_colors) +
  theme_bw() +
  facet_wrap(Alpha_div ~ ., scales = "free", nrow = 2, labeller = labeller(Alpha_div= Alpha_div_custom_labels)) +
  # RT category strip
  geom_segment(data = strip_labels_df,
               aes(x = x_start, xend = x_end, y = -0.5, yend = -0.5),
               inherit.aes = FALSE,
               size = 1.2) +
  geom_text(data = strip_labels_df,
            aes(x = label_pos, y = -0.8, label = RT_category),
            inherit.aes = FALSE,
            size = 5,
            fontface = "bold") +
  coord_cartesian(clip = "off")  + # allow drawing outside plot area
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22, face = "bold"),
        title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text = element_text(size = 18, face = "bold")) 

ggsave("./Plots/Alpha_div_boxplots_by_sampletype.png", plot=savePlot, height = 14, width = 24, units = "cm", limitsize = FALSE, dpi = 150)

# Add age grp to the above plot
ggplotObj <- ggplot(Merged_df_long %>% filter(SampleTypev2 != 'Other' ), aes(x = Alpha_div, y = Div_Value, fill = SampleTypev2, colour = AgeGroup))


savePlot <- ggplotObj + 
  geom_boxplot(position = position_dodge2(preserve = "single",padding = 0.3)) + # position = position_dodge(1.5)
  labs(title = "Alpha diversity for all samples (finer resolution)",
       fill = "Sample Type",
       x = NULL,
       y = "Diversity Value") + 
  scale_fill_manual(values = SampleType_custom_colors) +
  theme_bw() +
  facet_wrap(Alpha_div ~ ., scales = "free", nrow = 2, labeller = labeller(Alpha_div= Alpha_div_custom_labels)) +
  theme(axis.text.y = element_text(size=14), 
        axis.title= element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) 

ggsave("./Plots/Alpha_div_boxplots_by_sampletype_ageGrp.png", plot=savePlot, height = 20, width = 36, units = "cm", limitsize = FALSE, dpi = 150)


# Plot by age group
ggplotObj <- ggplot(Merged_df_long %>% filter(SampleTypev2 != 'Other' & Alpha_div %in% c("observed") & AgeGroup != 'NA'), aes(x = Alpha_div, y = Div_Value, color = AgeGroup))


savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity (Observed) for all samples",
       fill = "Sample Type",
       x = NULL,
       y = "Diversity Value") + 
  theme_bw() +
  facet_wrap(. ~ SampleTypev2, nrow = 2, scales = "free") +
  theme(axis.text.y = element_text(size=14), 
        axis.title= element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) 

ggsave("./Plots/Observed_diversity_boxplots_by_sampletype.png", plot=savePlot, height = 14, width = 30, units = "cm", limitsize = FALSE, dpi = 150)

# repeat above plot with selected sample types Buccal, sputum, saliva and BAL
# case 1
Disease_of_Interest <- c("Healthy", "Pneumonia", "Cystic Fibrosis")

ggplotObj <- Merged_df_long %>% 
  filter(SampleTypev2 %in% c("Oral_swab", "Nasal_Swab") & Alpha_div %in% c("diversity_shannon") & AgeGroup != 'NA' & Disease %in% Disease_of_Interest) %>% 
  ggplot(aes(x = Alpha_div, y = Div_Value, fill = Disease))


savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity (Shannon)",
       fill = "Health Status",
       x = NULL,
       y = "Diversity Value") + 
  theme_bw() +
  facet_wrap(AgeGroup ~ SampleTypev2, nrow = 1, scales = "free", labeller = labeller(.multi_line = F)  ) +
  theme(axis.text.y = element_text(size=14), 
        axis.title= element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) 


ggsave("./Plots/Observed_diversity_boxplots_DE_Pneumonia.png", plot=savePlot, height = 14, width = 16, units = "cm", limitsize = FALSE, dpi = 150)

# case 2
Disease_of_Interest <- c("Healthy", "Covid-19", "Cystic Fibrosis")

ggplotObj <- Merged_df_long %>% 
  filter(Alpha_div %in% c("Shannon") & AgeGroup != 'NA' & Disease %in% Disease_of_Interest & SampleTypev2 == "Sputum") %>% 
  ggplot(aes(x = Alpha_div, y = Div_Value, fill = Disease))


savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity (Shannon) in Sputum",
       fill = "Health Status",
       x = NULL,
       y = "Diversity Value") + 
  theme_bw() +
#  facet_wrap(AgeGroup ~ ., nrow = 1, scales = "free", labeller = labeller(.multi_line = F)  ) +
  theme(axis.text.y = element_text(size=14), 
        axis.title= element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) 


ggsave("./Plots/Observed_diversity_boxplots_DE_Sputum.png", plot=savePlot, height = 14, width = 20, units = "cm", limitsize = FALSE, dpi = 150)



# case 3 LRT
Disease_of_Interest <- c("Healthy", "Pneumonia")

ggplotObj <- Merged_df_long %>% 
  filter(RT_category %in% c("URT") & Alpha_div %in% c("diversity_shannon") & AgeGroup != 'NA' & Disease %in% Disease_of_Interest) %>% 
  ggplot(aes(x = Alpha_div, y = Div_Value, fill = Disease))


savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity (Shannon) in URT samples",
       fill = "Health Status",
       x = NULL,
       y = "Diversity Value") + 
  theme_bw() +
  #facet_wrap(AgeGroup ~ ., nrow = 1, scales = "free", labeller = labeller(.multi_line = F)  ) +
  theme(axis.text.y = element_text(size=14), 
        axis.title= element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) 


ggsave("./Plots/Observed_diversity_boxplots_DE_BAL.png", plot=savePlot, height = 14, width = 20, units = "cm", limitsize = FALSE, dpi = 150)

ggsave("./Plots/Observed_diversity_boxplots_DE_URT_withoutAge.png", plot=savePlot, height = 14, width = 20, units = "cm", limitsize = FALSE, dpi = 150)


# Rt category panel plot
ggplotObj <- ggplot(Merged_df_long %>% filter(SampleTypev2 != 'Other' & Alpha_div %in% c("observed") & AgeGroup != 'NA'), aes(x = Alpha_div, y = Div_Value, color = AgeGroup))



savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity (Observed) for all RT categories across age groups",
       fill = "Sample Type",
       x = NULL,
       y = "Diversity Value") + 
  theme_bw() +
  facet_wrap(. ~ RT_category, scales = "free", nrow = 1) +
  theme(axis.text.y = element_text(size=14), 
        axis.title= element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"))

ggsave("./Plots/Observed_diversity_boxplots_by_RT_cat.png", plot=savePlot, height = 14, width = 30, units = "cm", limitsize = FALSE, dpi = 150)

# Check corr between read count and alpha diversity
Merged_df_long$Total_read_count <- Merged_df_long$After_QC_R1 + Merged_df_long$After_QC_R2
Merged_df_long 

ggplotObj <- ggplot(Merged_df_long %>% filter(SampleTypev2 != 'Other' & Alpha_div %in% c("diversity_gini_simpson", "dominance_gini") & AgeGroup != 'NA')%>% group_by(SampleTypev2), aes(x = log10(Total_read_count), y = Div_Value, color = AgeGroup, linetype = Alpha_div))

savePlot <- ggplotObj + geom_smooth(method = 'lm', formula = y ~ poly(x, 2)) + 
  geom_point(alpha = 0.5, position = "jitter", size = 1.5)+
  facet_wrap(RT_category ~ AgeGroup , scales = "free_y", nrow = , labeller = label_wrap_gen(multi_line = F)) +
  labs(title = "Alpha diversity (y) v/s Read count log10 (x): y ~ poly(x, 2)",
       colour = "Alpha div",
       x = "Total read count (log10)",
       y = "Simpson diversity index") + 
  theme_bw() +
  xlim(4,9) +
  theme(axis.text.y = element_text(size=14), 
        axis.title= element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        axis.text.x = element_text(size = 12),
        #axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"))
  

ggsave("./Plots/Simpson_diversity_by_ReadCount_RT_cat.png", plot=savePlot, height = 20, width = 36, units = "cm", limitsize = FALSE, dpi = 150)








