# Script to do add apha diversity cols to metadata df
library(tidyverse)

# read Metadata table
Metadata_df <- read.table("../MiPoRT_Metadata_global_QC_passed_M4_profiled.tsv", header = T, sep = "\t")

# remove blank rows
Metadata_df <- Metadata_df[rowSums(is.na(Metadata_df)) < ncol(Metadata_df), ]

# List all diversity file names
Alpha_div_files <- c("../MetaPhlan4_results/Diversity_Merged/Metaphlan4_taxonomic_profiles_merged_gini.tsv", 
                     "../MetaPhlan4_results/Diversity_Merged/Metaphlan4_taxonomic_profiles_merged_richness.tsv",
                     "../MetaPhlan4_results/Diversity_Merged/Metaphlan4_taxonomic_profiles_merged_shannon.tsv",
                     "../MetaPhlan4_results/Diversity_Merged/Metaphlan4_taxonomic_profiles_merged_simpson.tsv")

Merged_df <- Metadata_df
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

# Remove samples with failed metaphaln diversity
Failed_Samples <- Merged_df %>% select("Gini.dominance":"Gini.simpson") %>% apply(MARGIN = 1, FUN = function(eachRow){
  sum(is.na(eachRow))
})

# subset_metadata_samples to remove metaphlan failed samples
Merged_df_sub <- Merged_df[-which(Failed_Samples > 0),]

# Update the same metadata table with alpha diversity values
write.table(Merged_df, "../MiPoRT_Metadata_global_QC_passed_M4_profiled.tsv", row.names = F, sep = '\t', quote = F)

# Reassign vars
Merged_df_trueCopy <- Merged_df_sub
Merged_df <- Merged_df_sub

# Let's keep nares for alpha div
#Merged_df <- Merged_df_sub %>% filter(SampleType != "Anterior_nares")


# Add a new sampletype col which only includes top 8 sampletypes
# Count frequency of each SampleType
top_sampletypes <- Merged_df %>%
  count(SampleType, sort = TRUE) %>%   # Count and sort
  slice_max(n, n = 9) %>%              # Get top 8 most frequent
  pull(SampleType)                     # Extract SampleType names


# Re-assign remaining sampletypes as "Other"
Merged_df <- Merged_df %>%
  mutate(SampleTypev2 = ifelse(as.character(SampleType) %in% top_sampletypes, as.character(SampleType), "Other"))

table(Merged_df$SampleTypev2, useNA = 'ifany')


# Plot diversity distribution
Merged_df_long <- Merged_df %>% pivot_longer(cols = c("Gini.dominance":"Gini.simpson"), names_to = "Alpha_div", values_to = "Div_Value")

# Add factors for SampleTypev2
SamplingSite_Factor <- c("Anterior_nares", "Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Sputum", "BAL", "Other")
SamplingSite_labels <- c("Anterior nares","Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL", "Other")

Merged_df_long$SampleTypev2 <- factor(Merged_df_long$SampleTypev2, 
                                 levels = SamplingSite_Factor,
                                 labels = SamplingSite_labels)  

levels(Merged_df_long$SampleTypev2)
# sanity check
table(Merged_df_long$SampleTypev2, useNA = 'ifany')

# Add factors for RT category
table(Merged_df_long$RT_category, useNA = 'ifany')
RT_Cat_Factor <- c("URT", "IRT", "LRT")
RT_Cat_Labels <- c("Upper RT", "Intermediate RT", "Lower RT")
Merged_df_long$RT_category <- factor(Merged_df_long$RT_category , 
                                 levels = RT_Cat_Factor,
                                 labels = RT_Cat_Labels)

# sanity check
table(Merged_df_long$RT_category, useNA = 'ifany')

# Add factors for health status
table(Merged_df_long$Healthy, useNA = 'ifany')
Health_stat_lvls <- c("TRUE", "FALSE", "Unknown")
Health_stat_labs <- c("Healthy (N=1480)", "Diseased (N=1685)", "Unknown (N=10)")
Merged_df_long$Healthy <- factor(Merged_df_long$Healthy, 
                             levels = Health_stat_lvls,
                             labels = Health_stat_labs)
# sanity check
table(Merged_df_long$Healthy, useNA = 'ifany')

# Define custom colors for top 9 + "Other"
SampleType_custom_colors <- c(
  "Sputum" = "#E41A1C", "BAL" = "#377EB8", "Tongue dorsum" = "#4DAF4A",
  "Supraglottal" = "#984EA3", "Oral swab" = "#FF7F00",
  "Nasal swab" = "#FFFF33", "Saliva" = "#A65628", "Buccal mucosa" = "#F781BF", "Anterior nares" = "#00CED1",
  "Other" = "#999999" )

# dim 3175*4 = 12700

# filter out Unknown health status samples
Merged_df_long <- Merged_df_long %>% filter(Healthy != "Unknown")
Merged_df_long$Healthy <- droplevels(Merged_df_long$Healthy)

# Plot 1. Alpha diversity density
ggplotObj <- ggplot(Merged_df_long, aes(Div_Value, fill = Healthy))

savePlot <- ggplotObj + 
  geom_density(kernel = "gaussian", alpha = 0.8) + 
  facet_wrap(Alpha_div ~ ., scales = "free", nrow = 2) +
  labs(title = "Alpha diversity distribution for health status", fill = "Health status") + 
  theme_bw() +
  theme(axis.text=element_text(size=22), 
        axis.title=element_text(size=22,face="bold"),
        title = element_text(size=24,face="bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text = element_text(size = 18, face = "bold"))  

# save plot
ggsave("./Plots/Alpha_div_dist_byHealth.png", 
       plot = savePlot,
       width = 35,              # Reduce dimensions
       height = 20,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE
)

# Plot 2 freq RT cat
# Custom labels for facets
Alpha_div_custom_labels <- c(
  "Gini.dominance" = "Gini_dominance",
  "Richness" = "Richness",
  "Shannon" = "Shannon",
  "Gini.simpson" = "Gini_simpson"
)

ggplotObj <- ggplot(Merged_df_long, aes(Div_Value, colour = RT_category))

savePlot <- ggplotObj + 
  geom_freqpoly(linewidth = 1.5) + 
  facet_wrap(~ Alpha_div, 
             nrow = 2, ncol = 2, 
             scales = "free", 
             labeller = labeller(Alpha_div = Alpha_div_custom_labels)) +
  labs(
    title = "Alpha diversity distribution for all samples (N=3175)",
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
    strip.text = element_text(size = 18, face = "bold")
  )

# save plot
ggsave("./Plots/Alpha_div_dist_by_RT.png", 
       plot = savePlot,
       width = 35,              # Reduce dimensions
       height = 20,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE
)

# Plot 3 freq with sample types
ggplotObj <- ggplot(Merged_df_long, aes(Div_Value, colour = SampleTypev2))

savePlot <- ggplotObj + geom_freqpoly(linewidth = 1.5) + 
  facet_wrap(Alpha_div ~ RT_category, scales = "free", nrow = 4, labeller = labeller(Alpha_div= Alpha_div_custom_labels,.multi_line = F)) +
  scale_color_manual(values = SampleType_custom_colors) +
  labs(title = "Alpha diversity distribution for all samples (N=3175)",
       colour = "Sample type (top 9)",
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
    strip.text = element_text(size = 18, face = "bold")
  )

# save plot
ggsave("./Plots/Alpha_div_dist_by_SampleTypev2.png",
       plot = savePlot,
       width = 50,          
       height = 30,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE)

# Plot 4: Box plot diversity by RT category
ggplotObj <- ggplot(Merged_df_long, aes(x = Alpha_div, y = Div_Value, fill = RT_category))

# save plot
savePlot <- ggplotObj + geom_boxplot() +
  labs(title = "Alpha diversity distribution for all samples (N=2968)",
       colour = "QC status",
       x = NULL,
       y = "Diversity Value") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22, face = "bold"),
        title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        strip.text = element_text(size = 18, face = "bold")) + 
  facet_wrap(Alpha_div ~ ., scales = "free", nrow = 2, labeller = labeller(Alpha_div= Alpha_div_custom_labels))

ggsave("./Plots/Alpha_div_boxplots_by_RT.png",
       plot = savePlot,
       width = 30,          
       height = 18,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE)

# Plot 5: Boxplot diversity by sampletypev2
table(Merged_df$SampleTypev2)

# SKIP Add factors Age groups
table(Merged_df$AgeGroup, useNA = 'ifany')
'      Adult       Child Older_adult Young_adult        <NA> 
        269        1074         432         802         391
'

AgeGrp_Factor<- c("Child", "Young_adult", "Adult", "Older_adult", "NA")
Merged_df_long$AgeGroup <- factor(Merged_df_long$AgeGroup, levels = AgeGrp_Factor)
levels(Merged_df_long$AgeGroup)

# plot by RT cat
ggplotObj <- ggplot(Merged_df_long %>% filter(SampleTypev2 != 'Other' ), aes(x = Alpha_div, y = Div_Value, fill = SampleTypev2))

savePlot <- ggplotObj + 
  geom_boxplot() +
  labs(title = "Alpha diversity for major sample types",
       fill = "Sample Type",
       x = NULL,
       y = "Diversity Value") + 
  scale_fill_manual(values = SampleType_custom_colors) +
  theme_bw() +
  facet_wrap(Alpha_div ~ ., scales = "free", nrow = 2, labeller = labeller(Alpha_div= Alpha_div_custom_labels)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 22, face = "bold"),
        title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        plot.margin = margin(t = 20, r = 20, b = 80, l = 20, unit = "pt"),
        strip.text = element_text(size = 18, face = "bold")) 

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








