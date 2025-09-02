# Script to do plot beta diversity statistics and boxplots within sample types
library(tidyverse)
library(introdataviz)
library(ggpubr)

# This matrix has samples from Anterior nares
Diversity_df <- read.table(
  "../04_calculate_diversity_postBC/MiPORT_All_RT_featureTable_sp_filt_1e-4_Relab_per_Min10_batchCorrected_samples_bray-curtis.tsv", 
  sep = '\t', header = T, check.names = F)

# Save colnames
SampleNames <- colnames(Diversity_df)

# Metadata file with sample profiles which passed metaphlan 
Metadata_df <- read.table(
  "../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", 
  header = T, sep = "\t")

  # sanity check
  glimpse(Metadata_df)

###########################################################################
# PART 1 Boxplot between beta div within sample types (with all samples)
###########################################################################
  
  # sanity check
  table(Metadata_df$SampleTypev2, useNA = 'ifany')
  
# subset metadata
Metadata_df_sub <- Metadata_df %>% 
  select(SampleID, RT_category, SampleTypev2, Healthy)

  # sanity check
  table(Metadata_df_sub$SampleTypev2, useNA = 'ifany')
  levels(Metadata_df_sub$SampleTypev2)

  # Make an ordered list for factors
  SamplingSite_Factor <- c("Anterior_nares", "Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Sputum", "BAL", "Other")

  # Format text labels
  SampleTypev2_labels <- c(
    "Anterior nares\n(n=208)",
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

  # For plotting define custom colors for top 9 + "Other"
  SampleType_custom_colors <- c(
    "Anterior nares\n(n=208)" = "#00BFC4",
    "Nasal swab\n(n=27)" = "#DAA520", 
    "Buccal mucosa\n(n=1027)" = "#F781BF",
    "Oral swab\n(n=246)" = "#FF7F00",
    "Saliva\n(n=315)" = "#A65628",
    "Tongue dorsum\n(n=418)" = "#4DAF4A",
    "Supraglottal\n(n=59)" = "#984EA3",
    "Sputum\n(n=193)" = "#E41A1C", 
    "BAL\n(n=578)" = "#377EB8",
    "Other\n(n=64)" = "#999999"  # Gray for "Other"
  )

  # Add factors
  Metadata_df_sub$SampleTypev2 <- factor(Metadata_df_sub$SampleTypev2,
                                       levels = SamplingSite_Factor,
                                       labels = SampleTypev2_labels)
  # sanity check
  levels(Metadata_df_sub$SampleTypev2)
  table(Metadata_df_sub$SampleTypev2)

  # convert diversity matrix into long format; Results in row names changed into col1 and col names into col2 values; with their values added into Distance Column
  Diversity_df_long <- as.data.frame(as.table(as.matrix(Diversity_df))) %>%
    rename(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
    filter(Sample1 != Sample2)  # Exclude diagonal elements
  
  # Add 1 cols for grouping information for all sample IDs
  Diversity_df_long <- Diversity_df_long %>%
    # Add sampletype for first column samples
    left_join(Metadata_df_sub, 
              by = c("Sample1" = "SampleID")) %>% 
    rename(Sample1_sampleType = SampleTypev2) %>%
    # Add sampletype for second column samples
    left_join(Metadata_df_sub%>% 
                select(-RT_category),
              by = c("Sample2" = "SampleID")) %>% 
    rename(Sample2_sampleType = SampleTypev2)
  
  # Filter and keep only inter-group comparisons
  Diversity_df_long <- Diversity_df_long %>%
    filter(Sample1_sampleType == Sample2_sampleType) # exclude inter sampletype comparisons
  
  # Remove reverse redundant groups
  Diversity_df_long <- Diversity_df_long %>% 
    filter(Sample1 < Sample2)

  # Plot 1: Boxplot beta diversity between groups
  ggPlotObj <- Diversity_df_long %>% group_by(Sample1_sampleType,
                                              RT_category) %>%
    ggplot(aes(x = Sample1_sampleType, y = Distance, 
               fill = Sample1_sampleType)) 
  
  # save plot
  savePlot <- ggPlotObj +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot
    labs(
      title = 'Beta Diversity Variation Among Individuals per Respiratory Sample Type',
      Fill = "Sample type",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    ) +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 
  
  # save png
  ggsave("./Batch_corrected_Plots/Boxplot_Beta-Diversity-Variation_all_SampleTypes.png",
         plot=savePlot, 
         height = 22, 
         width = 50, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/Boxplot_Beta-Diversity-Variation_all_SampleTypes.svg",
         plot=savePlot, 
         height = 22, 
         width = 50, 
         units = "cm", 
         limitsize = FALSE)
  
  # Plot 2. Same box plot without Anterior nares
  # Boxplot beta diversity between groups
  ggPlotObj <- Diversity_df_long %>% 
    filter(!Sample1_sampleType %in% c("Anterior nares\n(n=208)")) %>%
    ggplot(aes(x = Sample1_sampleType, y = Distance, 
               fill = Sample1_sampleType)) 
  
  # save plot
  savePlot <- ggPlotObj +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot
    labs(
      title = 'Beta Diversity Variation Among Individuals per Respiratory Sample Type',
      Fill = "Sample type",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    ) +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 
  
  # save png
  ggsave("./Batch_corrected_Plots/wo_AN_Boxplot_Beta-Diversity-Variation_all_SampleTypes.png",
         plot=savePlot, 
         height = 22, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/wo_AN_Boxplot_Beta-Diversity-Variation_all_SampleTypes.svg",
         plot=savePlot, 
         height = 22, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE)
  
  # statistics
  # Format text labels
  Comparisons<- list(
      #  c("Anterior nares\n(n=208)", "BAL\n(n=578)"),
      c("Nasal swab\n(n=27)", "BAL\n(n=578)"),
      c("Buccal mucosa\n(n=1027)", "BAL\n(n=578)"),
      c("Oral swab\n(n=246)", "BAL\n(n=578)"),
      c("Saliva\n(n=315)", "BAL\n(n=578)"),
      c("Tongue dorsum\n(n=418)", "BAL\n(n=578)"),
      c("Supraglottal\n(n=59)", "BAL\n(n=578)"),
      c("Sputum\n(n=193)", "BAL\n(n=578)"),
      c("Nasal swab\n(n=27)", "Sputum\n(n=193)"),
      c("Buccal mucosa\n(n=1027)", "Sputum\n(n=193)"),
      c("Oral swab\n(n=246)", "Sputum\n(n=193)"),
      c("Saliva\n(n=315)", "Sputum\n(n=193)"),
      c("Tongue dorsum\n(n=418)", "Sputum\n(n=193)"),
      c("Supraglottal\n(n=59)", "Sputum\n(n=193)")
    )

  # Compare groups pairwise
  pairwise_result <- pairwise.wilcox.test(Diversity_df_long$Distance, Diversity_df_long$Sample1_sampleType, p.adjust.method = "BH")
  print(pairwise_result)

  # save results
  write.table(data.frame(pairwise_result[["p.value"]], check.names = F),
            "Compare_betadiversity_between_sampleTypes.txt",
            row.names = T,
            sep = '\t')

  # Kruskal_result test for >3 group comparisons
  kruskal_result <- kruskal.test(Distance ~ Sample1_sampleType, data = Diversity_df_long)
  
  print(kruskal_result)
  print(paste("Overall group difference is significant:", kruskal_result$p.value < 0.005))
  
  # Boxplot with stats
  # Boxplot beta diversity between groups. Subset AN and Other sample types
  ggPlotObj <- Diversity_df_long %>% 
    filter(!Sample1_sampleType %in% c("Anterior nares\n(n=208)",
                                      "Other\n(n=64)")) %>%
    ggplot(aes(x = Sample1_sampleType, y = Distance, 
               fill = Sample1_sampleType))
  
  savePlot <- ggPlotObj +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot
    labs(
      title = 'Beta Diversity Variation Among Individuals per Respiratory Sample Type',
      subtitle = 'Wilcoxon comparisons vs BAL & Sputum; BH-adjusted p-values displayed (ns hidden).',
      Fill = "Sample type",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    )  +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) +
    stat_compare_means(comparisons = Comparisons,
                       method = "wilcox.test", # Statistical test
                       label = "p.format", # Use "p.signif" for significance stars, or "p.format" for exact p-values
                       p.adjust.method = "BH",
                       hide.ns = TRUE)

  # save png
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_all_SampleTypes.png",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_all_SampleTypes.svg",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE)
  
  # save one plot with stars
  savePlot <- ggPlotObj +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot
    labs(
      title = 'Comparison of beta diversity variation among individuals per respiratory sample type',
      subtitle = "Wilcoxon comparisons vs BAL & Sputum (BH-adjusted)\tSignificance: * p<0.05, ** p<0.01, *** p<0.001, **** p<1e-10",
      Fill = "Sample type",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    )  +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          plot.subtitle=element_text(size=16,face="italic"),
          plot.title = element_text(size=20,face="bold"),
          ) +
    stat_compare_means(comparisons = Comparisons,
                       method = "wilcox.test", # Statistical test
                       label = "p.signif", # Use "p.signif" for significance stars, or "p.format" for exact p-values
                       p.adjust.method = "BH",
                       hide.ns = TRUE)

  # save png
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_all_SampleTypes_withStars.png",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_all_SampleTypes_withStars.svg",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE)

###########################################################################
# PART 2 Box plot between beta div within sample types (with Health status)
###########################################################################

# Add Healthy state information to the distance box plots
  # Add metadata info and filter
  # convert diversity matrix into long format; Results in row names changed into col1 and col names into col2 values; with their values added into Distance Column
  Diversity_df_long <- as.data.frame(as.table(as.matrix(Diversity_df))) %>%
    rename(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
    filter(Sample1 != Sample2)  # Exclude diagonal elements
  
  # Add cols for grouping information for all sample IDs
  Diversity_df_long <- Diversity_df_long %>%
    # Add sampletype for first column samples
    left_join(Metadata_df_sub, 
              by = c("Sample1" = "SampleID")) %>% 
    rename(Sample1_sampleType = SampleTypev2,
           Sample1_Health = Healthy) %>%
    # Add sampletype for second column samples
    left_join(Metadata_df_sub%>% 
                select(-RT_category),
              by = c("Sample2" = "SampleID")) %>% 
    rename(Sample2_sampleType = SampleTypev2,
           Sample2_Health = Healthy) %>%
    # Filter for both Healthy or both diseased
    filter(Sample1_Health == Sample2_Health)
  
    # sanity check
  table(Diversity_df_long$Sample1_Health, 
        Diversity_df_long$Sample2_Health, 
        useNA = 'ifany')
  
  # Filter and keep only inter-group comparisons
  Diversity_df_long <- Diversity_df_long %>%
    filter(Sample1_sampleType == Sample2_sampleType) # exclude inter sampletype comparisons
  
  # Remove reverse redundant groups
  Diversity_df_long <- Diversity_df_long %>% 
    filter(Sample1 < Sample2)
  
  # sanity check
  table(Diversity_df_long$Sample1_Health, 
        Diversity_df_long$Sample2_Health, 
        useNA = 'ifany')
  
  # Add factors to Healthy
  levels(Diversity_df_long$Sample1_Health)
  
  Diversity_df_long$Sample1_Health <-
    factor(Diversity_df_long$Sample1_Health,
           levels = c('TRUE', 'FALSE', "Unknown"),
           labels = c('Healthy', 'Diseased', "Unknown"))
  
  # set colors for plots
  Grp_Healthy_custom_color <- c(
    "Healthy" = "#59D7D9",  # Sky blue
    "Diseased" = "#E6550D"  # A rich burnt orange
  )
  
  # remove sample types with Unknown health status & smaples present in only one group
  Diversity_df_long_sub <- Diversity_df_long %>% 
    filter(Sample1_Health != 'Unknown',
           !Sample1_sampleType %in% 
             c("Anterior nares\n(n=208)",
               "Tongue dorsum\n(n=418)",
               "Supraglottal\n(n=59)",
               "Other\n(n=64)"                                      ))
  
  # drop unused levels
  Diversity_df_long_sub$Sample1_Health <- 
      droplevels(Diversity_df_long_sub$Sample1_Health)
    # sanity check
    levels(Diversity_df_long_sub$Sample1_Health)
    table(Diversity_df_long_sub$Sample1_Health)

  Diversity_df_long_sub$Sample1_sampleType <-
    droplevels(Diversity_df_long_sub$Sample1_sampleType)
  # sanity check
  levels(Diversity_df_long_sub$Sample1_sampleType)
    table(Diversity_df_long$Sample1_sampleType)  

    # Statistics
    
    # Compare groups pairwise
    Compare_disease_result <- compare_means(
      Distance ~ Sample1_Health,
      Diversity_df_long_sub,
      group.by = "Sample1_sampleType", 
      p.adjust.method = "BH")
    
    Compare_disease_result %>% 
      dplyr::arrange(p.adj)
    
    # Make a plot df
    pval_df <- Compare_disease_result %>%
      dplyr::mutate(label = ifelse(p.adj < 1e-10, "****", 
                                   ifelse(p.adj < 0.001, "***", 
                                          ifelse(p.adj < 0.01, "**",
                                                 ifelse(p.adj < 0.05, "*", "ns"))))) %>%
      dplyr::rename(.y = Sample1_sampleType ) %>%
      dplyr::mutate(y.position = 
                      max(Diversity_df_long_sub$Distance, na.rm = TRUE) + 0.1)  # Adjust as needed
  
    # save statistics
    
    # save results
    write.table(pval_df,
                "Compare_betadiversity_between_Health_in_each_sampleTypes.txt",
                row.names = F,
                sep = '\t')
    
  # Plot 1: Boxplot beta diversity between groups
  ggPlotObj <- Diversity_df_long_sub %>% 
    ggplot(aes(x = Sample1_sampleType, y = Distance, 
               fill = Sample1_Health,
               color = Sample1_Health))

  savePlot <- ggPlotObj +
    introdataviz::geom_split_violin(alpha = .4, 
                                    trim = FALSE, 
                                    na.rm = TRUE) +
    geom_boxplot(width = .4, alpha = .5, fatten = NULL, show.legend = F)  +  # Boxplot
    labs(
      title = 'Comparison of beta diversity variation among individuals\nper respiratory sample type',
      subtitle = "Wilcoxon comparisons (BH-adjusted)",
      fill = "Health status",
      color = "Health status",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    )  +
    stat_summary(fun.data = "mean_cl_boot", 
                 geom = "pointrange", 
                 color = 'whitesmoke',
                 show.legend = F, 
                 position = position_dodge(.40)) +
    facet_wrap(~ Sample1_sampleType, scales = "free_x") + 
    theme_bw() +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          plot.subtitle=element_text(size=16,face="italic"),
          plot.title = element_text(size=20,face="bold"),
          plot.caption = element_text(size=14,face="bold"),
          strip.text = element_blank(), # rm facet wrap labels
          strip.background = element_blank(),
          ) +
  scale_color_manual(values = Grp_Healthy_custom_color) +
  scale_fill_manual(values = Grp_Healthy_custom_color) +
    # Compare Healthy vs Diseased within each sample type
    stat_compare_means(aes(label = after_stat(p.signif)),
                       label = "p.format",
                       method = "wilcox.test",
                       hide.ns = TRUE,
                       vjust = -5,
                       show.legend = FALSE)
  
' # did not work  
stat_pvalue_manual(pval_df,
  label = "label",
  y.position = "y.position",
                         xmin = "group1",
                         xmax = "group2",
                         tip.length = 0.01,
                         facet = ".y")'
  
  # save png
  ggsave("./Batch_corrected_Plots/wo_AN_Health_and_Stats_Boxplot_Beta-Diversity.png",
         plot=savePlot, 
         height = 30, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/wo_AN_Health_and_Stats_Boxplot_Beta-Diversity.svg",
         plot=savePlot, 
         height = 30, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)
  
  # save same plot with significance stars
  savePlot <- ggPlotObj +
    introdataviz::geom_split_violin(alpha = .4, 
                                    trim = FALSE, 
                                    na.rm = TRUE) +
    geom_boxplot(width = .4, alpha = .5, fatten = NULL, show.legend = F)  +  # Boxplot
    labs(
      title = 'Comparison of beta diversity variation among individuals per\nrespiratory sample type',
      subtitle = "Wilcoxon comparisons (BH-adjusted)",
      caption = "Significance: * p<0.05, ** p<0.01, *** p<0.001, **** p<1e-10",
      fill = "Health status",
      color = "Health status",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    )  +
    stat_summary(fun.data = "mean_cl_boot", 
                 geom = "pointrange", 
                 color = 'whitesmoke',
                 show.legend = F, 
                 position = position_dodge(.40)) +
    facet_wrap(~ Sample1_sampleType, scales = "free_x") + 
    theme_bw() +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          plot.subtitle=element_text(size=16,face="italic"),
          plot.title = element_text(size=20,face="bold"),
          strip.text = element_blank(), # rm facet wrap labels
          strip.background = element_blank(),
    ) +
    scale_color_manual(values = Grp_Healthy_custom_color) +
    scale_fill_manual(values = Grp_Healthy_custom_color) +
    # Compare Healthy vs Diseased within each sample type
    stat_compare_means(aes(label = after_stat(p.signif)),
                       label = "p.signif",
                       method = "wilcox.test",
                       hide.ns = TRUE,
                       vjust = -5,
                       show.legend = FALSE)
  # save png
  ggsave("./Batch_corrected_Plots/wo_AN_Health_and_Stats_Boxplot_Beta-Diversity_with_Stars.png",
         plot=savePlot, 
         height = 30, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/wo_AN_Health_and_Stats_Boxplot_Beta-Diversity_with_Stars.svg",
         plot=savePlot, 
         height = 30, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)

###########################################################################
# PART 3 Box plot between beta div within sample types (only Healthy)
###########################################################################
  
  # subset metadata for this work
  Metadata_df_sub <- Metadata_df %>% 
    select(SampleID, RT_category, SampleTypev2, Healthy) %>%
    filter(Healthy == 'TRUE',
           SampleTypev2 %in% SamplingSite_Factor) 
  
  # sample types to rm
  Samples_to_rm <- c("Anterior_nares", "Supraglottal", "Other")
  
  # Make an ordered list for remaining sampletypes for factors
  SamplingSite_Factor <- c("Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Sputum", "BAL")
  
  # Correct the sample sizes
  table(Metadata_df_sub$SampleTypev2, useNA = 'ifany')
  
  '
   BAL Buccal_mucosa    Nasal_Swab     Oral_swab        Saliva 
   50           527            18           171           210 
   Sputum Tongue_dorsum 
       17           418
  '
  # Format text labels
  SampleTypev2_labels <- c(
    #"Anterior nares\n(n=208)",
    "Nasal swab\n(n=18)",
    "Buccal mucosa\n(n=527)",
    "Oral swab\n(n=171)",
    "Saliva\n(n=210)",
    "Tongue dorsum\n(n=418)",
    #"Supraglottal\n(n=59)",
    "Sputum\n(n=17)",
    "BAL\n(n=50)"
    #"Other\n(n=64)"
  )
  levels(Metadata_df_sub$SampleTypev2)
  table(Metadata_df_sub$SampleTypev2)
  
  # subset Diversity df based on samples in metadata
  Diversity_df_sub <- Diversity_df[Metadata_df_sub$SampleID, # subset rows
                                   Metadata_df_sub$SampleID] # subset cols
  
    # sanity check
    dim(Diversity_df_sub)
  
  # convert diversity matrix into long format; Results in row names changed into col1 and col names into col2 values; with their values added into Distance Column
  Diversity_df_long <- as.data.frame(as.table(as.matrix(Diversity_df_sub))) %>%
    rename(Sample1 = Var1, Sample2 = Var2, Distance = Freq) %>%
    filter(Sample1 != Sample2)  # Exclude diagonal elements
  
  # sanity check. Should be less than this number of rows
  1411^2
  
  # Add cols for grouping information for all sample IDs
  Diversity_df_long <- Diversity_df_long %>%
    # Add sampletype for first column samples
    left_join(Metadata_df_sub, 
              by = c("Sample1" = "SampleID")) %>% 
    rename(Sample1_sampleType = SampleTypev2,
           Sample1_Health = Healthy) %>%
    # Add sampletype for second column samples
    left_join(Metadata_df_sub%>% 
                select(-RT_category),
              by = c("Sample2" = "SampleID")) %>% 
    rename(Sample2_sampleType = SampleTypev2,
           Sample2_Health = Healthy) 
  
  # dims should be 1990921 > 1989510 
  # sanity check
  table(Diversity_df_long$Sample1_Health, 
        Diversity_df_long$Sample2_Health, 
        useNA = 'ifany')
  
  # Filter and keep only intra-sampletype comparisons
  Diversity_df_long <- Diversity_df_long %>%
    filter(Sample1_sampleType == Sample2_sampleType) # exclude inter sampletype comparisons
    # dims are 1990921 > 1989510 > 527496
  
  # Remove reverse redundant groups
  Diversity_df_long <- Diversity_df_long %>% 
    filter(Sample1 < Sample2)
    # dims are 1990921 > 1989510 > 527496 > 263748
  
  # sanity check
  table(Diversity_df_long$Sample1_Health, 
        Diversity_df_long$Sample2_Health, 
        useNA = 'ifany')
  
  table(Diversity_df_long$Sample1_sampleType, 
        Diversity_df_long$Sample2_sampleType, 
        useNA = 'ifany')
  '
                     BAL Buccal_mucosa Nasal_Swab Oral_swab Saliva Sputum Tongue_dorsum
  BAL             1225             0          0         0      0      0             0
  Buccal_mucosa      0        138601          0         0      0      0             0
  Nasal_Swab         0             0        153         0      0      0             0
  Oral_swab          0             0          0     14535      0      0             0
  Saliva             0             0          0         0  21945      0             0
  Sputum             0             0          0         0      0    136             0
  Tongue_dorsum      0             0          0         0      0      0         87153
  
  '
  # Add factors to sample type
  table(Diversity_df_long$Sample1_sampleType)
  levels(Diversity_df_long$Sample1_sampleType)
  
  # check first
  SamplingSite_Factor
  SampleTypev2_labels
  
  # add factors now
  Diversity_df_long$Sample1_sampleType <- 
    factor(Diversity_df_long$Sample1_sampleType,
           levels = SamplingSite_Factor,
           labels = SampleTypev2_labels)
  
  # sanity check
  levels(Diversity_df_long$Sample1_sampleType)
  
  # now set custom colors for each sample type
  # For plotting define custom colors for top 9 + "Other"
  SampleType_custom_colors <- c(
    "Anterior nares\n(n=208)" = "#00BFC4",
    "Nasal swab\n(n=18)" = "#DAA520", 
    "Buccal mucosa\n(n=527)" = "#F781BF",
    "Oral swab\n(n=171)" = "#FF7F00",
    "Saliva\n(n=210)" = "#A65628",
    "Tongue dorsum\n(n=418)" = "#4DAF4A",
    "Supraglottal\n(n=59)" = "#984EA3",
    "Sputum\n(n=17)" = "#E41A1C", 
    "BAL\n(n=50)" = "#377EB8",
    "Other\n(n=64)" = "#999999"  # Gray for "Other"
  )
  
  # Add stats grouping to compare the bet div between BAL and sputum v/s others
  # statistics
  # Format text labels
  Comparisons<- list(
    #  c("Anterior nares\n(n=208)", "BAL\n(n=50)"),
    c("Nasal swab\n(n=18)", "BAL\n(n=50)"),
    c("Buccal mucosa\n(n=527)", "BAL\n(n=50)"),
    c("Oral swab\n(n=171)", "BAL\n(n=50)"),
    c("Saliva\n(n=210)", "BAL\n(n=50)"),
    c("Tongue dorsum\n(n=418)", "BAL\n(n=50)"),
    c("Sputum\n(n=17)", "BAL\n(n=50)"),
    c("Nasal swab\n(n=18)", "Sputum\n(n=17)"),
    c("Buccal mucosa\n(n=527)", "Sputum\n(n=17)"),
    c("Oral swab\n(n=171)", "Sputum\n(n=17)"),
    c("Saliva\n(n=210)", "Sputum\n(n=17)"),
    c("Tongue dorsum\n(n=418)", "Sputum\n(n=17)")
  )
  
  # Compare groups pairwise
  pairwise_result <- 
    pairwise.wilcox.test(Diversity_df_long$Distance, 
                         Diversity_df_long$Sample1_sampleType, 
                         p.adjust.method = "BH")
  
  print(pairwise_result)
  View(data.frame(pairwise_result[["p.value"]]))
  
  # save results
  write.table(data.frame(pairwise_result[["p.value"]], check.names = F),
              "Healthy_Compare_betadiversity_between_sampleTypes.txt",
              row.names = T,
              sep = '\t')
  
  # Kruskal_result test for >3 group comparisons
  kruskal_result <- 
    kruskal.test(Distance ~ Sample1_sampleType, 
                 data = Diversity_df_long)
  
  print(kruskal_result)
  print(paste("Overall group difference is significant:", kruskal_result$p.value < 0.005))
  
  # Boxplot beta diversity between groups.
  ggPlotObj <- Diversity_df_long %>%
    ggplot(aes(x = Sample1_sampleType, 
               y = Distance, 
               fill = Sample1_sampleType))
  
  savePlot <- ggPlotObj +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot
    labs(
      title = 'Beta diversity variation among healthy individuals per respiratory sample type',
      subtitle = 'Wilcoxon comparisons vs BAL & Sputum; BH-adjusted p-values displayed.',
      Fill = "Sample type",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    )  +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) +
    stat_compare_means(comparisons = Comparisons,
                       method = "wilcox.test", # Statistical test
                       label = "p.format", # Use "p.signif" for significance stars, or "p.format" for exact p-values
                       p.adjust.methods = "BH",
                       hide.ns = TRUE)
  
  # save png
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_Healthy_SampleTypes.png",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_Healthy_SampleTypes.svg",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE)

  # same plot but with stars
  savePlot <- ggPlotObj +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot
    labs(
      title = 'Beta diversity variation among healthy individuals per respiratory sample type',
      subtitle = 'Wilcoxon comparisons vs BAL & Sputum; BH-adjusted p-values displayed.',
      Fill = "Sample type",
      x = "Sample type (N)",
      y = "Bray-Curtis Distance"
    )  +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) +
    stat_compare_means(comparisons = Comparisons,
                       method = "wilcox.test", # Statistical test
                       label = "p.signif", # Use "p.signif" for significance stars, or "p.format" for exact p-values
                       p.adjust.methods = "BH",
                       hide.ns = FALSE)
  
  # save png
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_Healthy_SampleTypes_withStars.png",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("./Batch_corrected_Plots/wo_AN_Stats_Boxplot_Beta-Diversity-Variation_Healthy_SampleTypes_withStars.svg",
         plot=savePlot, 
         height = 35, 
         width = 45, 
         units = "cm", 
         limitsize = FALSE)
  