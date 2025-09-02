# Script to do PCoA on filtered and Metadata merged, BC similarity metrics
library(tidyverse)
library(vegan)

# This matrix has samples from Anterior nares
DistanceMetric_df <- read.table(
  "./MiPORT_All_RT_featureTable_sp_filt_1e-4_Relab_per_Min10_batchCorrected_samples_bray-curtis.tsv", 
  sep = '\t', header = T, check.names = F)

# Save colnames
SampleNames <- colnames(DistanceMetric_df)

# Metadata file with sample profiles which passed metaphlan 
Metadata_df <- read.table(
  "../../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", 
  header = T, sep = "\t")
  # sanity check
  glimpse(Metadata_df)
  
  # remove blank rows
  Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]
  
  # Rm QC failed samples if present
  table(Metadata_df$QC_Status_R1, Metadata_df$QC_Status_R2)
  # All 3135 samples retained.

# 
  table(Metadata_df$Healthy, useNA = 'ifany')
  # 3145 + 27 samples + metadata cols
  '
  FALSE    TRUE Unknown 
   1457    1668      10 
  '

# Ordination scaled by sample sizes. Less sample count should also be represented.
  # sanity check
  table(Metadata_df$SampleType, useNA = 'ifany')

  # Total sample size per SampleType
  sampletype_counts <- table(Metadata_df$SampleType, useNA = 'ifany')

# Calculate weights based on the reciprocal of sample size
weights <- case_when(
  Metadata_df$SampleType == "Anterior_nares" ~ 1 / sampletype_counts["Anterior_nares"],
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

# Perform PCoA on the dissimilarity matrix
bray_curtis_dist <- as.dist(DistanceMetric_df)

# Specify number of components
nComponents <- 3
pcoa_result <- wcmdscale(bray_curtis_dist, 
                         eig = TRUE, 
                         k = nComponents, 
                         w = weights)
  # sanity checks
  # Total number of axes with positive eigenvalues
  length(pcoa_result$eig[pcoa_result$eig > 0])  
  
  # Check if any eigenvalues are negative (non-Euclidean artifacts)
  summary(pcoa_result$eig)

# Calculate the variance explained by each PC
variance_explained <- pcoa_result$eig / sum(pcoa_result$eig) * 100  
  # Convert to percentages

# Add variance explained to PC names
PC_names <- c(paste0("PCoA-", 1:nComponents, 
                     " (", round(variance_explained[1:nComponents], 3), "%)"))

# Extract the coordinates for the PCs
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- PC_names

# {Outlier detection} Compute Z-scores for PC1 & PC2
  # Compute Z-scores for PC1 & PC2
    pcoa_coords$PC1_z <- c(scale(pcoa_coords$`PCoA-1 (14.041%)`))
    pcoa_coords$PC2_z <- c(scale(pcoa_coords$`PCoA-2 (6.359%)`))
    pcoa_coords$PC3_z <- c(scale(pcoa_coords$`PCoA-3 (5.91%)`))

  # Define an outlier threshold (e.g., |z| > 2.5 or 3)
  outlier_threshold <- 2.5
  pcoa_coords$outlier <- (abs(pcoa_coords$PC1_z) > outlier_threshold) | 
    (abs(pcoa_coords$PC2_z) > outlier_threshold) 
  # | (abs(pcoa_coords$PC3_z) > outlier_threshold)

# count outliers
  table(pcoa_coords$outlier, useNA = 'ifany')
    '
    FALSE  TRUE 
    3133     2'
  
# Extract outliers for labeling
  outliers <- pcoa_coords[pcoa_coords$outlier, ]

# left join data matrix and sample IDs
  Merged_df <- pcoa_coords %>% 
  rownames_to_column(var="SampleID") %>% 
  left_join(Metadata_df, by = join_by("SampleID" == "SampleID"))

# sanity check
  table(Merged_df$SampleType, useNA = "ifany")

# Add factors

  # 1. Add factors Sample type
  SamplingSite_Factor <- c("Anterior_nares","Nasal_Swab", "Nasopharyngeal_Aspirate", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Palatine_Tonsils","Throat", "Sputum", "BAL")
  SamplingSite_Labels <- c("Anterior nares","Nasal swab", "Nasopharyngeal aspirate", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Palatine tonsils","Throat", "Sputum", "BAL")
  
  # Without nares
  #SamplingSite_Factor <- c("Nasal_Swab", "Nasopharyngeal_Aspirate", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Palatine_Tonsils","Throat", "Sputum", "BAL")
  
  #SamplingSite_Labels <- c("Nasal swab", "Nasopharyngeal aspirate", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Palatine tonsils","Throat", "Sputum", "BAL")
  
  # sanity check
  length(unique(Merged_df$SampleType))

  # Add factors
  Merged_df$SampleType <- factor(Merged_df$SampleType, 
                                 levels = SamplingSite_Factor,
                                 labels = SamplingSite_Labels)
  
  # sanity check
  levels(Merged_df$SampleType)

  # Count frequency of each SampleType
  top_sampletypes <- Merged_df %>%
    count(SampleType, sort = TRUE) %>%   # Count and sort
    slice_max(n, n = 8) %>%              # Get top 8 most frequent
    pull(SampleType)                     # Extract SampleType names
  
  # drop unwanted levels
  top_sampletypes <- droplevels(top_sampletypes)
    levels(top_sampletypes)
    
    #Merged_df$SampleType <- droplevels(Merged_df$SampleType)
    #levels(Merged_df$SampleType)

  # Add a new sample type col which only includes top 8 sample types
  # Re-assign remaining sampletypes as "Other"
  Merged_df <- Merged_df %>%
    mutate(SampleTypev2 = ifelse(as.character(SampleType) %in% top_sampletypes, as.character(SampleType), "Other"))
  
  # sanity check
  table(Merged_df$SampleType, useNA = 'ifany')
  table(Merged_df$SampleTypev2, useNA = 'ifany')
  levels(Merged_df$SampleTypev2)
  
  # Add factors
  SampleTypev2_labels <- c(
    "Anterior nares (n=208)",
    #"Nasal swab (n=27)",
    "Buccal mucosa (n=1027)",
    "Oral swab (n=246)",
    "Saliva (n=315)",
    "Tongue dorsum (n=418)",
    "Supraglottal (n=59)",
    "Sputum (n=193)",
    "BAL (n=578)",
    "Other (n=91)"
  )
  
  # For plotting
  # Define custom colors for top 8 + "Other"
  SampleType_custom_colors <- c(
    "Anterior nares (n=208)" = "#00BFC4",
    "Nasal swab (n=27)" = "#DAA520", 
    "Buccal mucosa (n=1027)" = "#F781BF",
    "Oral swab (n=246)" = "#FF7F00",
    "Saliva (n=315)" = "#A65628",
    "Tongue dorsum (n=418)" = "#4DAF4A",
    "Supraglottal (n=59)" = "#984EA3",
    "Sputum (n=193)" = "#E41A1C", 
    "BAL (n=578)" = "#377EB8",
    "Other (n=91)" = "#999999"  # Gray for "Other"
  )
  
  Merged_df$SampleTypev2 <- 
    factor(Merged_df$SampleTypev2, 
           levels = c(levels(top_sampletypes), "Other"),
           labels = SampleTypev2_labels)
  table(Merged_df$SampleTypev2)
  levels(Merged_df$SampleTypev2)


# 2. Add factors RT category
  table(Merged_df$RT_category)
  
  RT_Cat_Factor <- c("URT", "IRT", "LRT")
  RT_Cat_Labels <- c("Upper RT (n=1838)", 
                     "Intermediate RT (n=719)", 
                     "Lower RT (n=578)")
  
  Merged_df$RT_category <- factor(Merged_df$RT_category, 
                                  levels = RT_Cat_Factor,
                                  labels = RT_Cat_Labels)
  
  # sanity check
  levels(Merged_df$RT_category)

# create a copy of Merged df
  Merged_df_trueCopy <- Merged_df
  # subset for plotting
  Merged_plot_df <- Merged_df %>%
    select(c("SampleID":"BioProject",
             RT_category,
             SampleTypev2
             ))
  
  # calculate centroids
  centroids_df <- Merged_plot_df %>%
    group_by(SampleTypev2) %>%
    summarise(
      x = mean(`PCoA-1 (14.041%)`, na.rm = TRUE),
      y = mean(`PCoA-2 (6.359%)`, na.rm = TRUE)
    )
  
  # centroids_df without Other sampletype
  centroids_df_sub <- centroids_df %>% filter(SampleTypev2 != 'Other (n=91)')
  centroids_df_sub$SampleTypev2 <- droplevels(centroids_df_sub$SampleTypev2)
  
  levels(centroids_df_sub$SampleTypev2)
  
  colnames(Merged_plot_df)
  '
   [1] "SampleID"         "PCoA-1 (14.041%)" "PCoA-2 (6.359%)" 
   [4] "PCoA-3 (5.91%)"   "PC1_z"            "PC2_z"           
   [7] "PC3_z"            "outlier"          "BioProject"      
  [10] "RT_category"      "SampleTypev2"
  '
# Plot 1
  # plot PC1 vs PC2
  ggplotObj <- Merged_plot_df %>% # filter(outlier == "FALSE") %>% 
    ggplot(aes(x=`PCoA-1 (14.041%)`, 
               y=`PCoA-2 (6.359%)`, 
               color=SampleTypev2, 
               fill = SampleTypev2 #, shape = outlier
               )
           )
  
# save plot
  savePlot <- ggplotObj + 
    geom_point(size = 1.5) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "PCoA based on Bray-Curtis\ndissimilarity (N=3135)", 
         colour = "Sample Type (N)") +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 
  
  savePlot # see plot & then save
  
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_sampletypev2.png", 
         plot=savePlot, 
         height = 18, 
         width = 26, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_sampletypev2.svg", 
         plot = savePlot, 
         height = 18, 
         width = 26, 
         units = "cm", 
         limitsize = FALSE)
  
  # save plot with ellipses without Other sampletype
  # plot PC1 vs PC2
  ggplotObj <- Merged_plot_df %>% 
    filter(SampleTypev2 != "Other (n=91)") %>% 
    ggplot(aes(x=`PCoA-1 (14.041%)`, 
               y=`PCoA-2 (6.359%)`, 
               color=SampleTypev2, 
               fill = SampleTypev2 #, shape = outlier
    )
    )
  
  savePlot <- ggplotObj + 
    geom_point(size = 1.5) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "PCoA based on Bray-Curtis dissimilarity (N=3135)", 
         subtitle = "Shaded ellipses represent 95% C.I. for each sample type", 
         colour = "Sample Type (N)") +
    stat_ellipse(aes(fill = SampleTypev2), geom = "polygon", alpha = 0.4, level = 0.95) +  # Confidence ellipses
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          subtitle = element_text(size=14, face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) +
    # To add centroids to each sample Type
    geom_point(data = centroids_df_sub, aes(x = x, y = y, color = SampleTypev2), 
                 shape = 4, size = 4, stroke = 2, show.legend = FALSE)
  #+ facet_wrap(SampleTypev2 ~ Healthy, nrow = 2)
  
  savePlot # see plot & then save
  
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_sampletypev2.png", 
         plot=savePlot, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_sampletypev2.svg", 
         plot = savePlot, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)
  
  #### plot PC2 vs PC3 ####
  ggplotObj <- Merged_plot_df %>% # filter(outlier == "FALSE") %>% 
    ggplot(aes(x=`PCoA-3 (5.91%)`, 
               y=`PCoA-2 (6.359%)`, 
               color=SampleTypev2, 
               fill = SampleTypev2)
#               shape = outlier
    )
  
  # save PC2 v/s PC3 plot without ellipses
  savePlot_1 <- ggplotObj + 
    geom_point(size = 1.5) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "PCoA based on Bray-Curtis \ndissimilarity (N=3135)", 
      colour = "Sample Type (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 
  
  savePlot_1 # see plot & then save
  
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_sampletypev2_02.png", 
         plot=savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_sampletypev2_02.svg", 
         plot = savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)
  
  # save PC2 v/s PC3 plot with ellipses
  savePlot_1 <- ggplotObj + 
    geom_point(size = 1.5) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "PCoA based on Bray-Curtis dissimilarity (N=3135)", 
         subtitle = "Shaded ellipses represent 95% C.I. for each sample type", 
         colour = "Sample Type (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) +
    stat_ellipse(aes(fill = SampleTypev2), geom = "polygon", alpha = 0.2, level = 0.95) +  # Confidence ellipses
    scale_fill_manual(values = SampleType_custom_colors)

  savePlot_1 # see plot & then save
  
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_sampletypev2_02.png", 
         plot=savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_sampletypev2_02.svg", 
         plot = savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)
  
  # plot PC1 vs PC3 without ellipses
  ggplotObj <- Merged_plot_df %>% # filter(outlier == "FALSE") %>% 
    ggplot(aes(x=`PCoA-1 (14.041%)`, 
               y=`PCoA-3 (5.91%)`, 
               color=SampleTypev2, 
               fill = SampleTypev2)
               #shape = outlier)
    )
  
  savePlot_2 <- ggplotObj + 
    geom_point(size = 1.5) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "PCoA based on Bray-Curtis\ndissimilarity (N=3135)", 
         colour = "Sample Type (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 

  
  savePlot_2 # see plot & then save
  
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_sampletypev2_03.png", 
         plot=savePlot_2, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_sampletypev2_03.svg", 
         plot = savePlot_2, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)

  # plot PC1 vs PC3 with ellipses
  ggplotObj <- Merged_plot_df %>% # filter(outlier == "FALSE") %>% 
    ggplot(aes(x=`PCoA-1 (14.041%)`, 
               y=`PCoA-3 (5.91%)`, 
               color=SampleTypev2, 
               fill = SampleTypev2)
           #shape = outlier)
    )
  
  savePlot_2 <- ggplotObj + 
    geom_point(size = 1.5) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "PCoA based on Bray-Curtis dissimilarity (N=3135)",
         subtitle = "Shaded ellipses represent 95% C.I. for each sample type", 
         colour = "Sample Type (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) +
    stat_ellipse(aes(fill = SampleTypev2), geom = "polygon", alpha = 0.2, level = 0.95) +  # Confidence ellipses
    scale_fill_manual(values = SampleType_custom_colors)
  
  
  savePlot_2 # see plot & then save
  
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_sampletypev2_03.png", 
         plot=savePlot_2, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_sampletypev2_03.svg", 
         plot = savePlot_2, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)

############ To save multiple plots ############

# Load library
  library(patchwork)

# Combine the plots into a panel (e.g., 1 row with 3 columns)
  combined_plot <- savePlot + savePlot_1 + savePlot_2 +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
  
  # watch and save
  combined_plot

  # Save combined plot
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_combined_panel.png",
       plot = combined_plot,
       height = 22, width = 54, units = "cm", dpi = 300)

  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_combined_panel.svg",
         plot = combined_plot,
         height = 22, width = 54, units = "cm")
  
############################################################
#### Repeat the 3 panel plot for RT category color now #####
############################################################

# Plot 2.1
# plot PC1 vs PC2
  # plot PC1 vs PC2
  ggplotObj <- Merged_plot_df %>% 
    filter(SampleTypev2 != "Other (n=91)") %>% 
    ggplot(aes(x=`PCoA-1 (14.041%)`, 
               y=`PCoA-2 (6.359%)`, 
               color=RT_category
    )
    )

    # save plot
    savePlot <- ggplotObj + 
      geom_point(size = 1.5) + 
      labs(title = "PCoA based on Bray-Curtis\ndissimilarity (N=3135)", 
           colour = "Respiratory category (N)") +
      theme_bw() + guides(fill="none") +
      theme(axis.text=element_text(size=20), 
            axis.title=element_text(size=20,face="bold"),
            title = element_text(size=20,face="bold"),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 18, face = "bold")) 
    # watch & save
    savePlot
    
    # save RT category plot
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_RT_CAT_01.png", 
       plot=savePlot, 
       height = 20, 
       width = 30, 
       units = "cm", 
       limitsize = FALSE, 
       dpi = 300)
    
    # save svg
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_RT_CAT_01.svg", 
       plot = savePlot, 
       height = 20, 
       width = 30, 
       units = "cm", 
       limitsize = FALSE)
    
    # save ordination with ellipses
    savePlot <- ggplotObj + 
      geom_point(size = 1.5) + 
      labs(title = "PCoA based on Bray-Curtis dissimilarity (N=3135)",
           subtitle = "Shaded ellipses represent 95% C.I. for each RT category", 
           colour = "Respiratory category (N)") +
      theme_bw() + guides(fill="none") +
      theme(axis.text=element_text(size=20), 
            axis.title=element_text(size=20,face="bold"),
            title = element_text(size=20,face="bold"),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 18, face = "bold")) +
      stat_ellipse(aes(fill = RT_category), geom = "polygon", alpha = 0.2, level = 0.95)
    # watch & save
    savePlot
    
    # save RT category plot
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_RT_CAT_01.png", 
           plot=savePlot, 
           height = 20, 
           width = 30, 
           units = "cm", 
           limitsize = FALSE, 
           dpi = 300)
    
    # save svg
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_RT_CAT_01.svg", 
           plot = savePlot, 
           height = 20, 
           width = 30, 
           units = "cm", 
           limitsize = FALSE)
    
# plot PC2 vs PC3
  # check PC names and update x & y
    PC_names

# plot obj
    ggplotObj <- Merged_plot_df %>% # filter(outlier == "FALSE") %>% 
    ggplot(aes(x=`PCoA-3 (5.91%)`, 
               y=`PCoA-2 (6.359%)`, 
               color=RT_category
               )
    )
  # save plot
  savePlot_1 <- ggplotObj + 
    geom_point(size = 1.5) + 
    labs(title = "PCoA based on Bray-Curtis\ndissimilarity (N=3135)", 
         colour = "Respiratory category (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 
  
  savePlot_1
  # save RT category plot
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_RT_CAT_02.png", 
         plot=savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_RT_CAT_02.svg", 
         plot = savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)
  
  # save PC2 vs PC3 ordination with ellipses
  savePlot_1 <- ggplotObj + 
    geom_point(size = 1.5) + 
    labs(title = "PCoA based on Bray-Curtis dissimilarity (N=3135)",
         subtitle = "Shaded ellipses represent 95% C.I. for each RT category", 
         colour = "Respiratory category (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) +
    stat_ellipse(aes(fill = RT_category), geom = "polygon", alpha = 0.2, level = 0.95)
  
  # watch and save
  savePlot_1
  
  # save RT category plot
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_RT_CAT_02.png", 
         plot=savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_RT_CAT_02.svg", 
         plot = savePlot_1, 
         height = 20, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE)
  
# plot PC1 vs PC3
  # check PC names and update x & y
  PC_names
  
  # create plot object
  ggplotObj <- Merged_plot_df %>% # filter(outlier == "FALSE") %>% 
  ggplot(aes(y=`PCoA-3 (5.91%)`, 
             x=`PCoA-1 (14.041%)`, 
             color=RT_category
             )
  )

  # save  PC1 vs PC3 ordination without ellipses
  savePlot_2 <- ggplotObj + 
    geom_point(size = 1.5) + 
    labs(title = "PCoA based on Bray-Curtis\ndissimilarity (N=3135)", 
         colour = "Respiratory category (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 
      
    # watch and save
    savePlot_2
    # save RT category plot
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_RT_CAT_03.png", 
           plot=savePlot_2, 
           height = 20, 
           width = 30, 
           units = "cm", 
           limitsize = FALSE, 
           dpi = 300)
    
    # save svg
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_RT_CAT_03.svg", 
           plot = savePlot_2, 
           height = 20, 
           width = 30, 
           units = "cm", 
           limitsize = FALSE)

    # save  PC1 vs PC3 ordination with ellipses
    savePlot_2 <- ggplotObj + 
      geom_point(size = 1.5) + 
      labs(title = "PCoA based on Bray-Curtis dissimilarity (N=3135)", 
           colour = "Respiratory category (N)") +
      theme_bw() + guides(fill="none") +
      theme(axis.text=element_text(size=20), 
            axis.title=element_text(size=20,face="bold"),
            title = element_text(size=20,face="bold"),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 18, face = "bold")) +
      stat_ellipse(aes(fill = RT_category), geom = "polygon", alpha = 0.2, level = 0.95)
    
    # watch and save
    savePlot_2
    # save RT category plot
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_RT_CAT_03.png", 
           plot=savePlot_2, 
           height = 20, 
           width = 30, 
           units = "cm", 
           limitsize = FALSE, 
           dpi = 300)
    
    # save svg
    ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_Ellipses_RT_CAT_03.svg", 
           plot = savePlot_2, 
           height = 20, 
           width = 30, 
           units = "cm", 
           limitsize = FALSE)
    
# To save multiple plots
  # Load library 
  #library(patchwork)

  # Combine the plots into a panel (e.g., 1 row with 3 columns)
  combined_plot <- savePlot + savePlot_1 + savePlot_2 +
    plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  
  # watch
  combined_plot
  
  # Save combined plot
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_combined_panel_RT_category.png",
         plot = combined_plot,
         height = 22, 
         width = 54, 
         units = "cm", 
         dpi = 300)

#### Side quest ####
  # plot PC1 v/s PC2 ordination with color sampletype and shapes with RT category

  # Save plot > plot PC1 vs PC2
  ggplotObj <- Merged_plot_df %>% # filter(outlier == "FALSE") %>% 
    ggplot(aes(x=`PCoA-1 (14.041%)`, 
               y=`PCoA-2 (6.359%)`, 
               color=SampleTypev2,
               shape = RT_category
    )
    )
  
  # save plot
  savePlot <- ggplotObj + 
    geom_point(size = 3) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "PCoA based on Bray-Curtis\ndissimilarity (N=3135)", 
         colour = "Sample Type (N)") +
    theme_bw() + guides(fill="none") +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) 
  
  savePlot
  # Save plot
  ggsave("Plots/PCoA_on_BC_diversity_BatchCorrected_SampleType_RT_category.png",
         plot = savePlot,
         height = 24, 
         width = 30, 
         units = "cm", 
         dpi = 300)
  
# 3D interactive plot
# Load required library
library(plotly)

  PC_names

  # Create 3D interactive ordination plot
  plot_ly(
  data = Merged_plot_df,
  x = ~`PCoA-1 (14.041%)`,
  y = ~`PCoA-2 (6.359%)`,
  z = ~`PCoA-3 (5.91%)`,
  color = ~SampleTypev2,
  colors = SampleType_custom_colors,
  text = ~paste("SampleID:", SampleID,
                "<br>RT Category:", RT_category,
                "<br>Sample Type:", SampleTypev2),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3)
) %>%
  layout(
    title = "3D Ordination Plot (Bray-Curtis PCoA)",
    scene = list(
      xaxis = list(title = "PCoA-1 (14.041%)"),
      yaxis = list(title = "PCoA-2 (6.359%)"),
      zaxis = list(title = "PCoA-3 (5.91%)")
    ),
    legend = list(font = list(size = 16))
  )

# save 3D plot
# Load required library
library(htmlwidgets)

# Create the sample type plot
ord_3d_plot <- plot_ly(
  data = Merged_plot_df,
  x = ~`PCoA-1 (14.041%)`,
  y = ~`PCoA-2 (6.359%)`,
  z = ~`PCoA-3 (5.91%)`,
  color = ~ SampleTypev2,
  colors = SampleType_custom_colors,
  text = ~paste("SampleID:", SampleID,
                "<br>RT Category:", RT_category,
                "<br>Sample Type:", SampleTypev2),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = list(
      text = "3D Ordination Plot (Bray-Curtis PCoA)",
      font = list(size = 22, family = "Arial", color = "black")
    ),
    margin = list(t = 100), 
    scene = list(
      xaxis = list(
        title = list(text = "PCoA-1 (14.041%)",
                     font = list(size = 22, family = "Arial", 
                                 color = "black")),
        tickfont = list(size = 18)
      ),
      yaxis = list(
        title = list(text = "PCoA-2 (6.359%)",
                     font = list(size = 22, family = "Arial", 
                                 color = "black")),
        tickfont = list(size = 18)
      ),
      zaxis = list(
        title = list(text = "PCoA-3 (5.91%)",
                     font = list(size = 22, family = "Arial", 
                                 color = "black")),
        tickfont = list(size = 18)
      )
    ),
    legend = list(
      title = list(text = "Sample Type", font = list(size = 18, 
                                                     face = "bold")),
      font = list(size = 22)
    )
  )

ord_3d_plot

# Save the plot as HTML
saveWidget(ord_3d_plot, file = "Batchcorrected_SampleType_3D_Ordination_Plot.html", selfcontained = TRUE)

# Plot with RT_category
# Create the plot
ord_3d_plot <- plot_ly(
  data = Merged_plot_df,
  x = ~`PCoA-1 (14.041%)`,
  y = ~`PCoA-2 (6.359%)`,
  z = ~`PCoA-3 (5.91%)`,
  color = ~RT_category,
  #colors = SampleType_custom_colors,
  text = ~paste("SampleID:", SampleID,
                "<br>RT Category:", RT_category,
                "<br>Sample Type:", SampleTypev2),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    title = list(
      text = "3D Ordination Plot (Bray-Curtis PCoA)",
      font = list(size = 22, family = "Arial", color = "black")
    ),
    margin = list(t = 100), 
    scene = list(
      xaxis = list(
        title = list(text = "PCoA-1 (14.041%)",
                     font = list(size = 22, family = "Arial", 
                                 color = "black")),
        tickfont = list(size = 18)
      ),
      yaxis = list(
        title = list(text = "PCoA-2 (6.359%)",
                     font = list(size = 22, family = "Arial", 
                                 color = "black")),
        tickfont = list(size = 18)
      ),
      zaxis = list(
        title = list(text = "PCoA-3 (5.91%)",
                     font = list(size = 22, family = "Arial", 
                                 color = "black")),
        tickfont = list(size = 18)
      )
    ),
    legend = list(
      title = list(text = "Sample Type", font = list(size = 18, 
                                                     face = "bold")),
      font = list(size = 22)
    )
  )

ord_3d_plot

# Save the plot as HTML
saveWidget(ord_3d_plot, file = "Batch_corrected_RT_category_3D_Ordination_Plot.html", selfcontained = TRUE)

########################################################################
 # plot without anterior nare samples
  
# plot 1.2 Outlier detection
  ggplotObj <- Merged_df %>% filter(Healthy != "Unknown" & SampleTypev2 %in% c("Sputum", "BAL")) %>% ggplot(aes(x=`Component_1 (12.83%)` , y=`Component_2 (8.91%)`, fill = Healthy, color=SampleTypev2, shape = outlier))
  
  savePlot <- ggplotObj + 
    geom_point(size = 1.5) + 
    scale_colour_manual(values = SampleType_custom_colors) +
    labs(title = "MDS based on BC dissimilarity (N=2863)", colour = "Sample Type") +
    stat_ellipse(aes(fill = Healthy), geom = "polygon", alpha = 0.2, level = 0.95) +  # Confidence ellipses
    theme_bw() +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold")) + facet_wrap(SampleTypev2 ~ ., nrow = 2)
  
  ########## rm 73 outliers
  Merged_df <- Merged_df %>% filter(outlier == "FALSE")
  
  # plot 2
  ggplotObj <- ggplot(Merged_df, aes(x=`Component_1 (12.83%)` , y=`Component_2 (8.91%)`, colour = RT_category))
  
  savePlot <- ggplotObj + 
    geom_point() + 
    labs(title = "MDS based on BC dissimilarity (N=2863)", colour = "RT category") +
    stat_ellipse(aes(fill = RT_category), geom = "polygon", alpha = 0.2, level = 0.95) +  # Confidence ellipses
    guides(fill="none") + # rm label for fill var
    theme_bw() + 
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold"))   
  
  ggsave("./BC_ordination_plots/MDS_filtered_BC_ellipse_RT_Category.png", plot=savePlot, height = 18, width = 24, units = "cm", limitsize = FALSE, dpi = 150)
  
  
  # plot 3
  # remove some missing values
  # Custom labels for facets
  Healthy_custom_labels <- c(
    "TRUE" = "Healthy",
    "FALSE" = "Diseased",
    "NA" = "Unknown"
  )
  
  which(is.na(Merged_df$Healthy))
  Merged_df <- Merged_df[-which(is.na(Merged_df$Healthy)),]
  
  Merged_df$Healthy <- factor(Merged_df$Healthy)
  levels(Merged_df$Healthy)

  
  # plot 3.1 RT category ellipse
  ggplotObj <- ggplot(Merged_df, aes(x=`Component_1 (12.83%)`, y=`Component_2 (8.91%)`, colour = RT_category))
  
  savePlot <- ggplotObj + 
    geom_point() + 
    labs(title = "MDS based on BC dissimilarity (N=2853)", colour = "Sample Type") +
    stat_ellipse(aes(fill = RT_category), geom = "polygon", alpha = 0.2, level = 0.95) +  # Confidence ellipses
    guides(fill="none") + # rm label for fill var
    theme_bw() +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 18, face = "bold")) +
    facet_grid(. ~ Healthy, labeller = labeller(Healthy=Healthy_custom_labels))
  
  ggsave("./BC_ordination_plots/PCoA_BC_filtered_Health_status_RT_Category_ellipse.png", plot=savePlot, height = 24, width = 36, units = "cm", limitsize = FALSE, dpi = 200)
  
  # plot 3.2 Age ellipse
  Merged_df$AgeGroup <- droplevels(Merged_df$AgeGroup)
  levels(Merged_df$AgeGroup)
  
  ggplotObj <- ggplot(Merged_df, aes(x=`Component_1 (12.83%)`, y=`Component_2 (8.91%)`, colour = AgeGroup))
  
  savePlot <- ggplotObj + 
    geom_point() + 
    labs(title = "MDS based on BC dissimilarity (N=2853)", colour = "Bioproject") +
    stat_ellipse(aes(fill = AgeGroup), geom = "polygon", alpha = 0.2, level = 0.95) +  # Confidence ellipses
    guides(fill="none") + # rm label for fill var
    theme_bw() +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 18, face = "bold")) +
    facet_grid(. ~ RT_category)
  
  ggsave("./BC_ordination_plots/MDS_filtered_BC_RT_cat_Dataset_preCorrection.png", plot=savePlot, height = 18, width = 30, units = "cm", limitsize = FALSE, dpi = 200)
  
  # plot 3.1 RT category ellipse
  
  SampleofInterest <- c("Buccal_mucosa", "Oral_swab", "Saliva", "Sputum", "BAL")
  
  ggplotObj <- Merged_df %>% filter(SampleTypev2 %in% SampleofInterest) %>% ggplot(aes(x=`Component_1 (12.83%)`, y=`Component_2 (8.91%)`, colour = SampleTypev2, shape= AgeGroup))
  
  savePlot <- ggplotObj + 
    geom_point() + 
    labs(title = "MDS based on BC dissimilarity (N=2853)", colour = "Sample type") +
    stat_ellipse(aes(fill = SampleTypev2), geom = "polygon", alpha = 0.2, level = 0.95) +  # Confidence ellipses
    guides(fill="none") + # rm label for fill var
    theme_bw() +
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 18, face = "bold")) +
    facet_grid(RT_category ~ Healthy, labeller = labeller(Healthy=Healthy_custom_labels))
  
  ggsave("./BC_ordination_plots/MDS_BC_filtered_Health_status_RT_Category__shape_ageGrp_Sampletype_ellipse.png", plot=savePlot, height = 34, width = 46, units = "cm", limitsize = FALSE, dpi = 200)
  # plot 4 
  Merged_df <- Merged_df_trueCopy
  Merged_df <- Merged_df[-which(is.na(Merged_df$Healthy)),]
  
  # subset for disease of interest
  table(Merged_df_trueCopy$Disease)
  
  # temp sol
  #Merged_df$Disease <- if_else(as.character(Merged_df$Disease) == "Healthy_genetically", "Healthy", as.character(Merged_df$Disease))
  
  Merged_df_sub <- Merged_df[Merged_df$Disease %in% c("Pneumonia", "Healthy", "Covid-19", "Cystic Fibrosis"),]
  
  # Add factors  
  Merged_df_sub$Disease <- factor(Merged_df_sub$Disease, levels = c("Healthy","Covid-19", "Cystic Fibrosis","Pneumonia"))
  
  table(Merged_df_sub$Disease)
  levels(Merged_df_sub$Disease)
  
  ggplotObj <- ggplot(Merged_df_sub, aes(x=`Component_1 (12.83%)`, y=`Component_2 (8.91%)`, colour = Disease, shape =Bioproject_plot))
  
savePlot <- ggplotObj + 
  geom_point() +
  #stat_ellipse(aes(fill = AgeGroup), geom = "polygon", alpha = 0.2, level = 0.95) +
  #guides(fill="none") + # rm label for fill var
  labs(title = "MDS based on BC dissimilarity (N=2347)", colour = "Disease", shape = "Dataset ID") + 
  theme_bw() +
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold"),
          strip.text.x = element_text(size = 16, face = "bold")) + # Column facet headers
    #strip.text.y = element_text(size = 18, face = "italic")) + # Row facet headers 
    facet_grid(. ~ RT_category)
  
  ggsave("./BC_ordination_plots/MDS_filtered_BC_Disease_RT_category_BatchEffect.png", plot=savePlot, height = 18, width = 32, units = "cm", limitsize = FALSE, dpi = 200)
  
  
# plot 5
  Merged_df <- Merged_df_trueCopy
  Merged_df$Total_Depth <- Merged_df$After_QC_R1 + Merged_df$After_QC_R2
  
  # Add a new column for Total_Depth (log10 scale)
  Merged_df$Depth_Category <- cut(
    log10(Merged_df$Total_Depth), 
    breaks = c(0, 5, 6, 7, 8, 9, 12),  # Customize these breaks for log10 fold differences
    labels = c("<10^5", "10^5-10^6", "10^6-10^7", "10^7-10^8", "10^8-10^9", ">=10^9"),
    include.lowest = TRUE
  )
  
  # Define custom colors for each category
    custom_colors <- c(
      "<10^5" = "#ffff33",  # Bright Red
      "10^5-10^6" = "#377eb8",  # Bright Blue
      "10^6-10^7" = "#4daf4a",  # Bright Green
      "10^7-10^8" = "#984ea3",  # Bright Purple
      "10^8-10^9" = "#ff7f00",  # Bright Orange
      ">=10^9" = "#e41a1c"  # Bright Yellow
    )
  
  # Create the ggplot object
  ggplotObj <- ggplot(Merged_df, aes(x = `Component_1 (12.83%)`, y = `Component_2 (8.91%)`, colour = Depth_Category))
  
  # Add points and color scale
  savePlot <- ggplotObj + 
    geom_point(size = 1.5) + 
    labs(title = "MDS based on BC dissimilarity (N=2853)", colour = "Read Depth") + 
    theme_bw() +
    theme(axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20, face = "bold"),
          title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold"),
          strip.text.x = element_text(size = 16, face = "bold")) +
    facet_grid(. ~ RT_category) +
    scale_colour_manual(values = custom_colors)
  
  # Save the plot
  ggsave("./BC_ordination_plots/MDS_filtered_BC_RT_category_readCount.png", plot=savePlot, height = 16, width = 28, units = "cm", limitsize = FALSE, dpi = 150)
  
  # Plot 6 sample types PCOAs
  # Loop plots part 1: Store all repeated measure project IDs > Group them for sampletypes > Iterate over them > PCA
  # Define custom colors for top 8 + "Other"
  SampleType_custom_colors <- c(
    "Sputum" = "#E41A1C", "BAL" = "#377EB8", "Tongue_dorsum" = "#4DAF4A",
    "Supraglottal" = "#984EA3", "Palatine_Tonsils" = "#FF7F00",
    "Anterior_nares" = "#FFFF33", "Saliva" = "#A65628", "Buccal_mucosa" = "#F781BF",
    "Throat" = "#999999"  # Gray for "Other"
  )
  
  # All repeated sample projects
  All_Sampletypes <- unique(Metadata_df$SampleType)
  # Big Projects only
  All_Sampletypes <- c("Buccal_mucosa", "Saliva", "Sputum", "BAL")
  
  SampleNames <- colnames(DistanceMetric_df)
  GGplot_Obj <- list()
  
  for (eachSampletype in All_Sampletypes) {
    #print(eachSampletype)
    # fetch samples IDs from specific project
    SampleIds_to_retain <- Metadata_df %>% filter(SampleType == eachSampletype) %>% pull(SampleID) %>% intersect(SampleNames)
    
    # Skip projects with no samples
    if (length(SampleIds_to_retain) == 0) {
      warning(paste("No samples found for project:", eachSampletype))
      next
    }
    
    # sanity check
    #print(length(SampleIds_to_retain))
    
    # subset BC distance matrix
    Subset_Distance_df <- DistanceMetric_df[SampleIds_to_retain, SampleIds_to_retain]
    print(paste("Subset df for", eachSampletype, "has following dims", nrow(Subset_Distance_df), ncol(Subset_Distance_df)))
    
    # Skip if not enough data for PCoA
    if (nrow(Subset_Distance_df) < 2) {
      warning(paste("Not enough samples for PCoA for sampletype:", eachSampletype))
      next
    }
    
    # Perform PCoA on the dissimilarity matrix
    bray_curtis_dist <- as.dist(Subset_Distance_df)
    
    nComponents <- 2
    pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = nComponents)
    
    # Calculate the variance explained by each PC
    variance_explained <- pcoa_result$eig / sum(pcoa_result$eig) * 100  # Convert to percentages
    
    # Add variance to component name
    PC_names <- c(paste0("Component_", 1:nComponents, 
                         " (", round(variance_explained[1:nComponents], 2), "%)"))
    
    # Extract the coordinates for the PCs
    pcoa_coords <- as.data.frame(pcoa_result$points)

    # left join data matrix and sample IDs
    Merged_df <- pcoa_coords %>% 
      rownames_to_column(var="SampleID") %>% 
      left_join(Metadata_df, by = join_by("SampleID" == "SampleID"))
    
    # Add factors RT category
    RT_Cat_Factor<- c("URT", "IRT", "LRT")
    Merged_df$RT_category <- factor(Merged_df$RT_category, levels = RT_Cat_Factor)
    
    # Add factor for Bio projects
    Merged_df$BioProject <- factor(Merged_df$BioProject)
    
    # Add factor for Bio projects
    Merged_df$AgeGroup <- factor(Merged_df$AgeGroup, levels = c("Child", "Young_adult", "Adult", "Older_adult", "NA"))
    
    # Plot colors check
    if (!all(levels(Merged_df$SampleType) %in% names(SampleType_custom_colors))) {
      stop("SampleType_custom_colors is missing some levels.")
    }
    
    # plot PCA
    TitlePlot <- paste("MDS for", eachSampletype, ": N=", nrow(Merged_df))
    ggplotObj <- ggplot(Merged_df, aes(x=V1, y=V2, colour = BioProject))
    
    savePlot <- ggplotObj + 
      geom_point(size = 1.5) + 
      #scale_colour_manual(values = SampleType_custom_colors) +
      stat_ellipse(aes(fill = BioProject), geom = "polygon", alpha = 0.2, level = 0.95) +
      guides(fill="none") + # rm label for fill var
      labs(title = TitlePlot, colour = "Project ID",
           x = PC_names[1], 
           y = PC_names[2]) + 
      theme_bw() +
      theme(axis.text=element_text(size=18), 
            axis.title=element_text(size=16,face="bold"),
            title = element_text(size=18,face="bold"),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 18, face = "bold"))
    
    # Save the ggplot object to the list
    GGplot_Obj[[eachSampletype]] <- savePlot
    
  }
  
  # Plot
  # Load required library
  library(patchwork)
  
  # Combine all plots in GGplot_Obj using +
  combined_plot <- wrap_plots(GGplot_Obj, 
                              nrow = 2) +  # Set number of columns
    plot_annotation(title = "MDS Plots across all projects by sample types",
                    theme = theme(plot.title = element_text(size = 18, face = "bold")))
  
  # Display the arranged plot
  print(combined_plot)
  
  ggsave("Repeated_Measures/Combined_samplesite_allProjects_bySampletype_Ellipse_Project_MDS.png", plot=combined_plot, height = 24, width = 50, units = "cm", limitsize = FALSE, dpi = 150)
  
  
# Plot 7: Loop plot: For all project IDs > Group them for Rt category > Iterate over them > PCA
  # Define custom colors for top 8 + "Other"
  SampleType_custom_colors <- c(
    "Sputum" = "#E41A1C", "BAL" = "#377EB8", "Tongue_dorsum" = "#4DAF4A",
    "Supraglottal" = "#984EA3", "Palatine_Tonsils" = "#FF7F00",
    "Anterior_nares" = "#FFFF33", "Saliva" = "#A65628", "Buccal_mucosa" = "#F781BF",
    "Throat" = "#999999"  # Gray for "Other"
  )
  
  # All repeated sample projects
  All_RTCats <- unique(Metadata_df$RT_category)
  # Big Projects only
  
  
  SampleNames <- colnames(DistanceMetric_df)
  GGplot_Obj <- list()
  
  for (eachRT_category in All_RTCats) {
    #print(eachRT_category)
    # fetch samples IDs from specific project
    SampleIds_to_retain <- Metadata_df %>% filter(RT_category == eachRT_category) %>% pull(SampleID) %>% intersect(SampleNames)
    
    # Skip projects with no samples
    if (length(SampleIds_to_retain) == 0) {
      warning(paste("No samples found for project:", eachRT_category))
      next
    }
    
    # sanity check
    #print(length(SampleIds_to_retain))
    
    # subset BC distance matrix
    Subset_Distance_df <- DistanceMetric_df[SampleIds_to_retain, SampleIds_to_retain]
    print(paste("Subset df for", eachRT_category, "has following dims", nrow(Subset_Distance_df), ncol(Subset_Distance_df)))
    
    # Skip if not enough data for PCoA
    if (nrow(Subset_Distance_df) < 2) {
      warning(paste("Not enough samples for PCoA for RT category:", eachRT_category))
      next
    }
    
    # Perform PCoA on the dissimilarity matrix
    bray_curtis_dist <- as.dist(Subset_Distance_df)
    
    nComponents <- 2
    pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = nComponents)
    
    # Extract the coordinates for the PCs
    pcoa_coords <- as.data.frame(pcoa_result$points)
    
    PC_names <- c(paste0("Component_", 1:nComponents))
    colnames(pcoa_coords) <- PC_names
    
    # left join data matrix and sample IDs
    Merged_df <- pcoa_coords %>% 
      rownames_to_column(var="SampleID") %>% 
      left_join(Metadata_df, by = join_by("SampleID" == "SampleID"))
    
    # Add factors RT category
    RT_Cat_Factor<- c("URT", "IRT", "LRT")
    Merged_df$RT_category <- factor(Merged_df$RT_category, levels = RT_Cat_Factor)
    
    # Add factor for Bio projects
    Merged_df$BioProject <- factor(Merged_df$BioProject)
    
    # Add factor for Bio projects
    Merged_df$AgeGroup <- factor(Merged_df$AgeGroup, levels = c("Child", "Young_adult", "Adult", "Older_adult", "NA"))
    
    # Plot colors check
    if (!all(levels(Merged_df$SampleType) %in% names(SampleType_custom_colors))) {
      stop("SampleType_custom_colors is missing some levels.")
    }
    
    # plot PCA
    TitlePlot <- paste("MDS for", eachRT_category, ": N=", nrow(Merged_df))
    ggplotObj <- ggplot(Merged_df[Merged_df$Healthy %in% c("TRUE", "FALSE"),], aes(x=Component_1, y=Component_2, colour = BioProject))
    
    savePlot <- ggplotObj + 
      geom_point(size = 1.5) + 
      #scale_colour_manual(values = SampleType_custom_colors) +
      labs(title = TitlePlot, colour = "Project ID") + 
      theme_bw() +
      theme(axis.text=element_text(size=18), 
            axis.title=element_text(size=16,face="bold"),
            title = element_text(size=18,face="bold"),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 18, face = "bold")) + facet_wrap(Healthy ~ .)
    
    # Save the ggplot object to the list
    GGplot_Obj[[eachRT_category]] <- savePlot
  }
  
  # Plot
  # Load required library
  #library(patchwork)
  
  # Combine all plots in GGplot_Obj using +
  combined_plot <- wrap_plots(GGplot_Obj, 
                              nrow = 2) +  # Set number of columns
    plot_annotation(title = "MDS Plots across all projects by respiratory tract category",
                    theme = theme(plot.title = element_text(size = 18, face = "bold")))
  
  # Display the arranged plot
  print(combined_plot)
  
  ggsave("BC_ordination_plots/Combined_Rtcat_allProjects_Healthy_MDS.png", plot=combined_plot, height = 24, width = 40, units = "cm", limitsize = FALSE, dpi = 150)
  