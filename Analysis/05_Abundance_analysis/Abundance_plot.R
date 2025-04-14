# Script to create a PA matrix from abundance profile and plot major group profiles
library(stringr)
library(tidyverse)

# Read feature table & metadata table
Metaphlan_profiles <- read.table("../MetaPhlan4_results/Metaphlan4_species_merged_filtered_batchCorrected.tsv", header = T, sep = "\t")
Metadata_Df <- read.table("../MiPoRT_Metadata_global_QC_passed_M4_profiled.tsv", header = T, sep = "\t")

# remove blank rows
Metadata_Df <- Metadata_Df[rowSums(is.na(Metadata_Df)) < ncol(Metadata_Df), ]
#Metadata_Df <- Metadata_Df %>% filter(SampleType != 'Anterior_nares')

# SKIP
# Add a new sampletype col
    # Count frequency of each SampleType
    top_sampletypes <- Metadata_Df %>%
        count(SampleType, sort = TRUE) %>%   # Count and sort
        slice_max(n, n = 9) %>%              # Get top 8 most frequent
        pull(SampleType)                     # Extract SampleType names
    
    # Re-assign remaining sampletypes as "Other"
    Metadata_Df <- Metadata_Df %>%
        mutate(SampleTypev2 = ifelse(as.character(SampleType) %in% top_sampletypes, as.character(SampleType), "Other"))
    
    write.table(Metadata_Df,"../MiPoRT_Metadata_global_QC_passed_M4_profiled.tsv", row.names = F, sep = "\t", quote = F)

# Add factors for SampleTypev2
SamplingSite_Factor <- c("Anterior_nares", "Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Sputum", "BAL", "Other")
SamplingSite_labels <- c("Anterior nares","Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL", "Other")

Metadata_Df$SampleTypev2 <- factor(Metadata_Df$SampleTypev2, 
                                   levels = SamplingSite_Factor,
                                   labels = SamplingSite_labels)  

levels(Metadata_Df$SampleTypev2)
# sanity check
table(Metadata_Df$SampleTypev2, useNA = 'ifany')

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
Health_stat_labs <- c("Healthy (N=1480)", "Diseased (N=1685)", "Unknown (N=10)")
Metadata_Df$Healthy <- factor(Metadata_Df$Healthy, 
                              levels = Health_stat_lvls,
                              labels = Health_stat_labs)
# sanity check
table(Metadata_Df$Healthy, useNA = 'ifany')

# Define custom colors for top 9 + "Other"
SampleType_custom_colors <- c(
    "Sputum" = "#E41A1C", "BAL" = "#377EB8", "Tongue dorsum" = "#4DAF4A",
    "Supraglottal" = "#984EA3", "Oral swab" = "#FF7F00",
    "Nasal swab" = "#FFFF33", "Saliva" = "#A65628", "Buccal mucosa" = "#F781BF", "Anterior nares" = "#00CED1",
    "Other" = "#999999" )

# Add factors for AgeGrp
AgeGrp_Factor<- c("Child", "Young_adult", "Adult", "Older_adult")
Metadata_Df$AgeGroup <- factor(Metadata_Df$AgeGroup, levels = AgeGrp_Factor)

levels(Metadata_Df$AgeGroup)

# Get a list of M4 profile samples to retain
retainSamples <- intersect((colnames(Metaphlan_profiles)), Metadata_Df$SampleID)

# subset metadata with this
Metadata_Df <- Metadata_Df %>% filter(SampleID %in% retainSamples)
Metaphlan_profiles <- Metaphlan_profiles %>% select(all_of(retainSamples))

# reorder rows to match cols of M4
Metadata_Df <- Metadata_Df[match(colnames(Metaphlan_profiles), Metadata_Df$SampleID), ]

# sanity check
identical(colnames(Metaphlan_profiles), Metadata_Df$SampleID)  # Should return TRUE

# Get number of samples
nSamples <- (dim(Metaphlan_profiles)[2])
nTaxa <- dim(Metaphlan_profiles)[1]

print(paste0("Check if these are the number of samples you have: ",nSamples))
print(paste0("Check if these are the number of taxa you have: ",nTaxa))

# Count frequency of non-zero species abundances
Zero_Status_Metaphlan_profiles <- Metaphlan_profiles == 0
# Convert logical TRUE (absence) and FALSE (presence) to numeric 0 and 1
Presence_Absence_Df <- as.data.frame(!Zero_Status_Metaphlan_profiles) * 1

# 1. Calculate global prevalence of each Taxa
Freq_Taxa <- data.frame(
    Species = rownames(Metaphlan_profiles),
    Global_Prevalence = round(rowMeans(Presence_Absence_Df) * 100, 3),  # Calculate % of non-zero values
    SamplePresence_count = rowSums(Presence_Absence_Df)
)

# Remove taxa with <10 sample presence
Freq_Taxa_min10_Samples <- Freq_Taxa[Freq_Taxa$SamplePresence_count > 10,]
# 2384 taxa > filter min 10 samples > 780 taxa

nTaxa <- dim(Freq_Taxa_min10_Samples)[1]

write.table(Presence_Absence_Df %>% filter(row.names(Presence_Absence_Df) %in% Freq_Taxa_min10_Samples$Species), "../MetaPhlan4_results/MiPORT_filtered_min10_samples_PA_data.txt", sep = "\t", quote = F, row.names = T)

# plot dist
ggObj <- ggplot(Freq_Taxa_min10_Samples, aes(Global_Prevalence))

# Histogram plot for prevalence
savePlot <- ggObj + geom_histogram(bins = 25) + theme_light(base_size = 14) + labs(title="Global prevalence dist of taxa (n=2935)", x= "Prevalence %", y= "Sample Count")
# save plot
ggsave(file="./Plots/Histogram_prevalence_taxa_min10.png", plot=savePlot, width = 20, height = 20, units = "cm")

# 2. Calculate RT_category wise abundance of each Taxa
# First we join metadata with Metaphlan_profiles by matching sample names (Sample names must match in both tables)
# Assuming `Profile_Df` column names are sample names and metadata has a column with sample names.

long_PA_df <- Metaphlan_profiles %>% rownames_to_column(var = "Taxonomy") %>% filter(Taxonomy %in% Freq_Taxa_min10_Samples$Species)  %>%
    pivot_longer(cols = c(2:(nSamples+1)), names_to = "SampleName", values_to = "Abundance") %>%
    left_join(Metadata_Df, by = c("SampleName" = "SampleID")) 

# sanity check 
table(long_PA_df$BioProject)/nTaxa
table(Metadata_Df$BioProject)

str(long_PA_df[,2:34])

# sanity check
levels(long_PA_df$RT_category)
table(long_PA_df$RT_category)/nTaxa

# calculate prevalence within RT cat
for(eachRTcat in c("Upper RT", "Intermediate RT", "Lower RT")){
    # create a subset for each sampleType
    Profile_subset <- long_PA_df %>% filter(RT_category == eachRTcat) %>% group_by(Taxonomy) 
    
    nSamples_subset <- nrow(Profile_subset)/nTaxa
    print(paste("Calculating abundance for samples from ", eachRTcat, "n = ", nSamples_subset))
    
    # calculate abundance in percentage of each Taxa
    Freq_Taxa_subset <- Profile_subset %>% 
        group_by(Taxonomy) %>%  
        summarize(Abundance = sum(Abundance)*100)
    
    ColName <- paste("Abundance", eachRTcat, sep = "_")
    colnames(Freq_Taxa_subset) <- c("Species", ColName)
    
    # Add abundance to existing df of prevalence
    Freq_Taxa_min10_Samples <- left_join(Freq_Taxa_min10_Samples, Freq_Taxa_subset, by = c("Species" = "Species"))
    #print(Freq_Taxa_subset)
}

# Order by Abundance in descending order
Freq_Taxa_min10_Samples <- Freq_Taxa_min10_Samples[order(Freq_Taxa_min10_Samples$`Abundance_Lower RT`, decreasing = TRUE), ]

# calculate prevalence within sampletypev2
SampleType <- unique(long_PA_df$SampleTypev2)[-10]

for(eachSampleType in SampleType){
    # create a subset for each sampleType
    Profile_subset <- long_PA_df %>% filter(SampleTypev2 == eachSampleType) %>% group_by(Taxonomy) 
    
    nSamples_subset <- nrow(Profile_subset)/nTaxa
    print(paste("Calculating abundance for samples from ", eachSampleType, "n = ", nSamples_subset))
    
    # calculate total abundance of each Taxa
    Freq_Taxa_subset <- Profile_subset %>% 
        group_by(Taxonomy) %>%  
        summarize(Prevalence = sum(Abundance)*100)
    
    ColName <- paste("Abundance", eachSampleType, sep = "_")
    colnames(Freq_Taxa_subset) <- c("Species", ColName)
    
    # Add prevalence to existing df of prevalence
    Freq_Taxa_min10_Samples <- left_join(Freq_Taxa_min10_Samples, Freq_Taxa_subset, by = c("Species" = "Species"))
    #print(Freq_Taxa_subset)
}

# rm Abundance_Anterior nares col

Freq_Taxa_min10_Samples <- select(Freq_Taxa_min10_Samples, -"Abundance_Anterior nares")

colnames(Freq_Taxa_min10_Samples)

# Melt for plotting
Freq_Taxa_min10_Samples_long_RT <- Freq_Taxa_min10_Samples %>% 
    select(all_of(c("Species", "Abundance_Upper RT", "Abundance_Intermediate RT", "Abundance_Lower RT"))) %>%
    pivot_longer(cols = c(2:4), names_to = "RT_category", values_to = "Abundance_RT_category") %>% 
    ungroup() %>% group_by(RT_category)

# export major table
write.table(Freq_Taxa_min10_Samples, "Abundance_Taxa_filtered_min10_samples_All_MajorGroups.txt", sep = '\t', row.names = F, quote = F)

# Melt
Freq_Taxa_min10_Samples_long <-  Freq_Taxa_min10_Samples %>% 
    select(all_of(c("Species", colnames(Freq_Taxa_min10_Samples)[7:14]))) %>%
    pivot_longer(cols = c(2:9), names_to = "SampleType", values_to = "Abundance_SampleType")

# Use Freq_Taxa_min10_Samples_long & Freq_Taxa_min10_Samples_long_RT to plot abundance now

# Plot 1: Plot top taxa prevalence in LRT
glimpse(Freq_Taxa_min10_Samples_long_RT)

levels(Freq_Taxa_min10_Samples_long_RT$RT_category)

# Add factors
Freq_Taxa_min10_Samples_long_RT$RT_category <- 
    factor(Freq_Taxa_min10_Samples_long_RT$RT_category, 
    levels = c("Abundance_Upper RT", "Abundance_Intermediate RT",  "Abundance_Lower RT" ))

# check
levels(Freq_Taxa_min10_Samples_long_RT$RT_category)
table(Freq_Taxa_min10_Samples_long_RT$RT_category, useNA = 'ifany')

# plot distribution
ggObj <- ggplot(Freq_Taxa_min10_Samples_long_RT, aes(log10(Abundance_RT_category), fill = RT_category))
savePlot <- ggObj + 
    geom_density(alpha = 0.6) + 
    labs(title="Total abundance of 780 species across RT", x= "Abundance log10", y= "Density") + 
    theme_bw() +
    theme(axis.text=element_text(size=22), 
          axis.title=element_text(size=22,face="bold"),
          title = element_text(size=24,face="bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold")) 

# save plot
ggsave("./Plots/Abundance_log10_dist_by_RT.png", 
       plot = savePlot,
       width = 30,              # Reduce dimensions
       height = 15,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE
)

# plot heatmap
# sort species with LRT prevalence

Species_Order <- Freq_Taxa_min10_Samples[order(Freq_Taxa_min10_Samples$`Abundance_Upper RT`, decreasing = F),] %>% pull(Species)

Freq_Taxa_min10_Samples_long_RT$Species <- factor(Freq_Taxa_min10_Samples_long_RT$Species,
                                                  levels = Species_Order,
                                                  ordered = T
)

ggObj <- ggplot(Freq_Taxa_min10_Samples_long_RT, aes(x=RT_category, y=Species, fill = log10(Abundance_RT_category)))

savePlot <- ggObj + geom_tile() + theme_light(base_size = 14) + 
    labs(x="RT category", title = "Total abundance for species across RT category", y= "Species (n=780)") +
    theme_bw() +
    theme(axis.text=element_text(size=22), 
          axis.title=element_text(size=22,face="bold"),
          title = element_text(size=24,face="bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold"),
          axis.text.y=element_blank(), #remove y axis labels
          axis.ticks.y=element_blank()  #remove y axis ticks
    ) + 
    scale_fill_gradientn(
        colours = c("#000000", "#1f78b4", "#ffffcc", "#feb24c", "#f03b20"),
        values = scales::rescale(c(-2, 0, 2, 3, 4)),
        limits = c(-2, 4),
        na.value = "grey90"
    ) 

# save plot
ggsave(file="./Plots/RT_category_Abundance_heatmap_Min10_samples.png", 
       plot = savePlot,
       width = 45,              # Reduce dimensions
       height = 30,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE
)

# Subset based on species of interest; Check file 

# select species based on LRT abundance
Taxa_of_Interest <- Freq_Taxa_min10_Samples[order(Freq_Taxa_min10_Samples$`Abundance_Lower RT`, decreasing = T),] %>% 
    top_n(25, wt = `Abundance_Lower RT`) %>% pull(Species)

Freq_Taxa_min10_Samples_long_RT <- Freq_Taxa_min10_Samples_long_RT %>% 
    mutate(SpeciesOfInterest = Species %in% Taxa_of_Interest) %>%
    filter(SpeciesOfInterest == TRUE)

Freq_Taxa_min10_Samples_long_RT_sub <- Freq_Taxa_min10_Samples_long_RT %>%
    filter(SpeciesOfInterest == "TRUE")

# Get order of sp
Species_Order <- Freq_Taxa_min10_Samples %>%
    filter(Species %in% Freq_Taxa_min10_Samples_long_RT_sub$Species) %>%
    arrange("Abundance_Upper RT")%>% 
    pull(Species) 

    # or select for URT
    Species_Order <- Freq_Taxa_min10_Samples %>%
        filter(Species %in% Freq_Taxa_min10_Samples_long_RT_sub$Species) %>%
        arrange("Abundance_Upper RT")%>% 
        pull(Species) 

    # Factor main df with this sp order
    Freq_Taxa_min10_Samples_long_RT$Species <- factor(Freq_Taxa_min10_Samples_long_RT$Species,
                                                          levels = Species_Order,
                                                          ordered = T
        )

# Factor main df with this sp order
Freq_Taxa_min10_Samples_long_RT_sub$Species <- factor(Freq_Taxa_min10_Samples_long_RT_sub$Species, 
                                                      levels = Species_Order
                                                      ,ordered = T)
# plot
ggObj <- ggplot(Freq_Taxa_min10_Samples_long_RT_sub, aes(x=RT_category, y=Species, fill = log10(Abundance_RT_category) ))

savePlot <- ggObj + geom_tile() + theme_light(base_size = 14) + 
    labs(x="RT category", title = "Top 25 abundant species in LRT", y= "Species") +
    theme_bw() +
    theme(axis.text=element_text(size=22), 
          axis.title=element_text(size=22,face="bold"),
          title = element_text(size=24,face="bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold"),
          #axis.text.y=element_blank(), #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
    ) +  
    scale_fill_gradientn(
        colours = c("#000000", "#1f78b4", "#ffffcc", "#feb24c", "#f03b20"),
        values = scales::rescale(c(-2, 0, 2, 3, 4)),
        limits = c(-2, 4),
        na.value = "grey90"
    ) 

# save plot
ggsave(file="./Plots/RT_category_Abundance_Top25_LRT_abundant_taxa.png",
       plot = savePlot,
       width = 60,              # Reduce dimensions
       height = 30,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE
)


# For the same list of taxa plot their abundance in sampletypes

str(Freq_Taxa_min10_Samples_long)

# add factors to sample types of interest
SamplingSite_labels <- c("Nasal swab", "Buccal mucosa", "Oral swab", "Saliva", "Tongue dorsum", "Supraglottal", "Sputum", "BAL")

Freq_Taxa_min10_Samples_long$SampleType <- factor(Freq_Taxa_min10_Samples_long$SampleType, 
levels = unique(Freq_Taxa_min10_Samples_long$SampleType),
labels = SamplingSite_labels,
ordered = T)

# sanity check
levels(Freq_Taxa_min10_Samples_long$SampleType)

# Add column for species of interest
Freq_Taxa_min10_Samples_long

# Add a column for species of interest to long df for plotting
Freq_Taxa_min10_Samples_long <- Freq_Taxa_min10_Samples_long %>% 
    mutate(SpeciesOfInterest = Species %in% Taxa_of_Interest) %>%
    filter(SpeciesOfInterest == TRUE)

Freq_Taxa_min10_Samples_long$SpeciesOfInterest <- factor(Freq_Taxa_min10_Samples_long$SpeciesOfInterest)

Freq_Taxa_min10_Samples_long_sub <- Freq_Taxa_min10_Samples_long %>%
    filter(SpeciesOfInterest == "TRUE")


# Get order of sp
Species_Order <- Freq_Taxa_min10_Samples %>%
    filter(Species %in% Freq_Taxa_min10_Samples_long_sub$Species) %>%
    arrange(`Abundance_Lower RT`)%>% 
    pull(Species) 

# Factor main df with this sp order
Freq_Taxa_min10_Samples_long_sub$Species <- factor(Freq_Taxa_min10_Samples_long_sub$Species, 
levels = Species_Order,
ordered = T)
# plot
ggObj <- ggplot(Freq_Taxa_min10_Samples_long_sub, aes(x=SampleType, y=Species, fill = log10(Abundance_SampleType) ))

savePlot <- ggObj + geom_tile() + theme_light(base_size = 14) + 
    labs(x="Sample Type", title = "Top 25 abundant species in BAL", y= "Species") +
    theme_bw() +
    theme(axis.text=element_text(size=22), 
          axis.title=element_text(size=22,face="bold"),
          title = element_text(size=24,face="bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold"),
          #axis.text.y=element_blank(), #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
    ) 

# save plot
ggsave(file="./Plots/SampleType_Top25_BALF_Abundant_taxa.png", 
              plot = savePlot,
              width = 100,              # Reduce dimensions
              height = 50,
              units = "cm",
              dpi = 300,               # High resolution
              limitsize = FALSE
       )

# Try boxplot
# plot
ggObj <- ggplot(Freq_Taxa_min10_Samples_long, aes(x=SampleType, y =log10(Abundance_SampleType), fill = SampleType ))

savePlot <- ggObj + 
    geom_boxplot(alpha = 0.5) +
    geom_point() +
    labs(x="Sample Type", title = "Abundance v/s Sampletype", y= "Top 25 Species in BAL") +
    scale_fill_manual(values = SampleType_custom_colors) +
    theme_bw() +
    theme(axis.text=element_text(size=22),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title=element_text(size=22,face="bold"),
          title = element_text(size=24,face="bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20, face = "bold"),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
          strip.text = element_text(size = 18, face = "bold"),
          #axis.text.y=element_blank(), #remove y axis labels
          #axis.ticks.y=element_blank()  #remove y axis ticks
    ) 

# save plot
ggsave(file="./Plots/Boxplot_SampleType_Top25_BALF_Abundant_taxa.png", 
       plot = savePlot,
       width = 40,              # Reduce dimensions
       height = 20,
       units = "cm",
       dpi = 300,               # High resolution
       limitsize = FALSE
)

# Boxplot for top species in URT










