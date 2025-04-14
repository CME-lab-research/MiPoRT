# Script to create a PA matrix from abundance profile and plot major group profiles
library(stringr)
library(tidyverse)

# Read feature table & metadata table
Metaphlan_profiles <- read.table("../MetaPhlan4_results/Metaphlan4_species_merged_filtered_batchCorrected.tsv", header = T, sep = "\t")
Metadata_Df <- read.table("../MiPORT_Metadata_downstream_filtered_M4_passed.tsv", header = T, sep = "\t")

# remove blank rows
Metadata_Df <- Metadata_Df[rowSums(is.na(Metadata_Df)) < ncol(Metadata_Df), ]
Metadata_Df <- Metadata_Df %>% filter(SampleType != 'Anterior_nares')

SamplingSiteFactor <- c("Nasal_Swab", "Nasopharyngeal_Aspirate", "Buccal_mucosa", "Oral_swab", "Saliva", "Tongue_dorsum", "Supraglottal", "Palatine_Tonsils","Throat", "Sputum", "BAL")
Metadata_Df$SampleType <- factor(Metadata_Df$SampleType, levels = SamplingSiteFactor)

table(Metadata_Df$SampleType)

# Count frequency of each SampleType
top_sampletypes <- Metadata_Df %>%
  count(SampleType, sort = TRUE) %>%   # Count and sort
  slice_max(n, n = 8) %>%              # Get top 8 most frequent
  pull(SampleType)                     # Extract SampleType names

top_sampletypes <- droplevels(top_sampletypes)
Merged_df$SampleType <- droplevels(Merged_df$SampleType)

levels(Merged_df$SampleType)

# Add a new sampletype col which only includes top 8 sampletypes
# Re-assign remaining sampletypes as "Other"
Metadata_Df <- Metadata_Df %>%
  mutate(SampleTypev2 = ifelse(as.character(SampleType) %in% top_sampletypes, as.character(SampleType), "Other"))

# Add factors
Metadata_Df$SampleTypev2 <- factor(Metadata_Df$SampleTypev2, levels = c(levels(top_sampletypes), "Other"))  

write.table(Metadata_Df,"../MiPORT_Metadata_downstream_filtered_M4_passed_v2.tsv", row.names = F, sep = "\t", quote = F)


levels(Metadata_Df$SampleTypev2)

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
  Global_Prevalence = round(rowMeans(Presence_Absence_Df) * 100, 3)  # Calculate % of non-zero values
)

# Remove taxa with <2 sample presence
Freq_Taxa_min2_Samples <- Freq_Taxa[Freq_Taxa$Global_Prevalence > (1/nSamples)*100,]
# 2384 taxa > filter min 2 samples > 1497 taxa

nTaxa <- dim(Freq_Taxa_min2_Samples)[1]

write.table(Presence_Absence_Df %>% filter(row.names(Presence_Absence_Df) %in% Freq_Taxa_min2_Samples$Species), "../MetaPhlan4_results/MiPORT_filtered_PA_data.txt", sep = "\t", quote = F, row.names = T)

  # plot dist
  ggObj <- ggplot(Freq_Taxa_min2_Samples, aes(Global_Prevalence))
  
  # Histogram plot for prevalence
  savePlot <- ggObj + geom_histogram(bins = 25) + theme_light(base_size = 14) + labs(title="Global prevalence dist of taxa (n=2935)", x= "Prevalence %", y= "Sample Count")
  # save plot
  ggsave(file="./Plots/Histogram_prevalence_taxa_min2.png", plot=savePlot, width = 20, height = 20, units = "cm")

# 2. Calculate RT_category wise prevalence of each Taxa
# First we join metadata with Metaphlan_profiles by matching sample names (Sample names must match in both tables)
# Assuming `Profile_Df` column names are sample names and metadata has a column with sample names.
  
long_PA_df <- Presence_Absence_Df %>% rownames_to_column(var = "Taxonomy") %>% filter(Taxonomy %in% Freq_Taxa_min2_Samples$Species)  %>%
  pivot_longer(cols = c(2:(nSamples+1)), names_to = "SampleName", values_to = "Abundance") %>%
  left_join(Metadata_Df, by = c("SampleName" = "SampleID")) 

# sanity check 
table(long_PA_df$BioProject)/nTaxa
table(Metadata_Df$BioProject)

str(long_PA_df[,2:34])
# 1. Add factors RT cat
long_PA_df$RT_category <- factor(long_PA_df$RT_category, levels = c("URT", "IRT", "LRT"))
# sanity check
levels(long_PA_df$RT_category)
table(long_PA_df$RT_category)/nTaxa

# calculate prevalence within RT cat
for(eachRTcat in c("URT", "IRT", "LRT")){
  # create a subset for each sampleType
  Profile_subset <- long_PA_df %>% filter(RT_category == eachRTcat) %>% group_by(Taxonomy) 
  
  nSamples_subset <- table(Profile_subset$RT_category)/nTaxa
  print(paste("Calculating prevalence for samples from ", eachRTcat, "n = ", nSamples_subset))
  
  # calculate prevalence in percentage of each Taxa
  Freq_Taxa_subset <- Profile_subset %>% 
    group_by(Taxonomy) %>%  
    summarize(Prevalence = sum(Abundance > 0)*100 / n())
  
  ColName <- paste("Prevalence", eachRTcat, sep = "_")
  colnames(Freq_Taxa_subset) <- c("Species", ColName)
  
  # Add prevalence to existing df of prevalence
  Freq_Taxa_min2_Samples <- left_join(Freq_Taxa_min2_Samples, Freq_Taxa_subset, by = c("Species" = "Species"))
  #print(Freq_Taxa_subset)
}

# Order by prevalence in descending order
Freq_Taxa_min2_Samples <- Freq_Taxa_min2_Samples[order(Freq_Taxa_min2_Samples$Prevalence_LRT, decreasing = TRUE), ]

# calculate prevalence within sampletypev2
SampleType <- unique(levels(long_PA_df$SampleTypev2))[-9]

for(eachSampleType in SampleType){
  # create a subset for each sampleType
  Profile_subset <- long_PA_df %>% filter(SampleTypev2 == eachSampleType) %>% group_by(Taxonomy) 
  
  nSamples_subset <- table(Profile_subset$SampleTypev2)/nTaxa
  print(paste("Calculating prevalence for samples from ", eachSampleType, "n = ", nSamples_subset))
  
  # calculate prevalence in percentage of each Taxa
  Freq_Taxa_subset <- Profile_subset %>% 
    group_by(Taxonomy) %>%  
    summarize(Prevalence = sum(Abundance > 0)*100 / n())
  
  ColName <- paste("Prevalence", eachSampleType, sep = "_")
  colnames(Freq_Taxa_subset) <- c("Species", ColName)
  
  # Add prevalence to existing df of prevalence
  Freq_Taxa_min2_Samples <- left_join(Freq_Taxa_min2_Samples, Freq_Taxa_subset, by = c("Species" = "Species"))
  #print(Freq_Taxa_subset)
}

# Melt for plotting
Freq_Taxa_min2_Samples_long_RT <- Freq_Taxa_min2_Samples %>% 
  select(all_of(c("Species", "Prevalence_URT", "Prevalence_IRT", "Prevalence_LRT"))) %>%
  pivot_longer(cols = c(2:4), names_to = "RT_category", values_to = "Prevalence_RT_category") %>% 
  ungroup() %>% group_by(RT_category)

# export major table
write.table(Freq_Taxa_min2_Samples, "Frequency_filtered_Taxa_All_MajorGroups.txt", sep = '\t', row.names = F, quote = F)

# Melt
Freq_Taxa_min2_Samples_long <-  Freq_Taxa_min2_Samples %>% 
  select(all_of(c("Species", colnames(Freq_Taxa_min2_Samples)[6:13]))) %>%
    pivot_longer(cols = c(2:9), names_to = "SampleType", values_to = "Prevalence_SampleType")

# Use Freq_Taxa_min2_Samples_long & Freq_Taxa_min2_Samples_long_RT to plot prevalence now

# Plot 1: Plot top taxa prevalence in LRT
str(Freq_Taxa_min2_Samples_long_RT)

# Add factors
Freq_Taxa_min2_Samples_long_RT$RT_category <- factor(Freq_Taxa_min2_Samples_long_RT$RT_category, levels = c("Prevalence_URT", "Prevalence_IRT", "Prevalence_LRT") )

# check
levels(Freq_Taxa_min2_Samples_long_RT$RT_category)

# plot distribution
ggObj <- ggplot(Freq_Taxa_min2_Samples_long_RT, aes(log2(Prevalence_RT_category), fill = RT_category))
savePlot <- ggObj + geom_density(alpha = 0.6) + theme_light(base_size = 14) + labs(title="Prevalence distribution by RT", x= "Prevalence log2", y= "Density")


# plot heatmap
# sort species with LRT prevalence

Species_Order <- Freq_Taxa_min2_Samples[order(Freq_Taxa_min2_Samples$Prevalence_URT, decreasing = F),] %>% pull(Species)

Freq_Taxa_min2_Samples_long_RT$Species <- factor(Freq_Taxa_min2_Samples_long_RT$Species,
                                                 levels = Species_Order,
                                                 ordered = T
                                                 )

ggObj <- ggplot(Freq_Taxa_min2_Samples_long_RT, aes(x=RT_category, y=Species, fill = Prevalence_RT_category ))

savePlot <- ggObj + geom_tile() + theme_light(base_size = 14) + 
  labs(x="Prevalence in RT category", title = "RT category", y= "Species (n=1497)") +
  theme_bw() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text.y=element_blank(), #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
        )

# save plot
ggsave(file="./Plots/RT_category_Prevalence_Min2_samples.png", plot=savePlot, width = 25, height = 12, units = "cm")

# Subset based on species of interest; Check file 
Taxa_of_Interest <- Freq_Taxa_min2_Samples[order(Freq_Taxa_min2_Samples$Prevalence_LRT, decreasing = T),] %>% 
  top_n(25, wt = Prevalence_LRT) %>% pull(Species)

Freq_Taxa_min2_Samples_long_RT <- Freq_Taxa_min2_Samples_long_RT %>% 
  mutate(SpeciesOfInterest = Species %in% Taxa_of_Interest) %>%
  filter(SpeciesOfInterest == TRUE)

Freq_Taxa_min2_Samples_long_RT_sub <- Freq_Taxa_min2_Samples_long_RT %>%
  filter(SpeciesOfInterest == "TRUE")


Species_Order <- Freq_Taxa_min2_Samples %>%
  filter(Species %in% Freq_Taxa_min2_Samples_long_RT_sub$Species) %>%
  arrange(Prevalence_URT)%>% 
  pull(Species) 

Freq_Taxa_min2_Samples_long_RT$Species <- factor(Freq_Taxa_min2_Samples_long_RT$Species,
                                                 levels = Species_Order,
                                                 ordered = T
)


# Get order of sp
Species_Order <- Freq_Taxa_min2_Samples %>%
  filter(Species %in% Freq_Taxa_min2_Samples_long_RT_sub$Species) %>%
  arrange(Prevalence_URT)%>% 
  pull(Species) 

# Factor main df with this sp order
Freq_Taxa_min2_Samples_long_RT_sub$Species <- factor(Freq_Taxa_min2_Samples_long_RT_sub$Species, 
                                                     levels = Species_Order
                                                     ,ordered = T)
# plot
ggObj <- ggplot(Freq_Taxa_min2_Samples_long_RT_sub, aes(x=RT_category, y=Species, fill = Prevalence_RT_category ))

savePlot <- ggObj + geom_tile() + theme_light(base_size = 14) + 
  labs(x="RT category", title = "Top 25 prevalent species in LRT", y= "Species") +
  theme_bw() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        #axis.text.y=element_blank(), #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
  )

# save plot
ggsave(file="./Plots/RT_category_Prevalence_Top25_LRT_prevalent_taxa.png", plot=savePlot, width = 30, height = 15, units = "cm")


# For the same list of taxa plot their prevalence in sampletypes

str(Freq_Taxa_min2_Samples_long)

# add factors to sample types of interest
Freq_Taxa_min2_Samples_long$SampleType <- factor(Freq_Taxa_min2_Samples_long$SampleType, 
                                                 levels = unique(Freq_Taxa_min2_Samples_long$SampleType),
                                                 ordered = T)

# sanity check
levels(Freq_Taxa_min2_Samples_long$SampleType)

# Add column for species of interest
Freq_Taxa_min2_Samples_long


Freq_Taxa_min2_Samples_long <- Freq_Taxa_min2_Samples_long %>% 
  mutate(SpeciesOfInterest = Species %in% Taxa_of_Interest) %>%
  filter(SpeciesOfInterest == TRUE)

Freq_Taxa_min2_Samples_long$SpeciesOfInterest <- factor(Freq_Taxa_min2_Samples_long$SpeciesOfInterest)

Freq_Taxa_min2_Samples_long_sub <- Freq_Taxa_min2_Samples_long %>%
  filter(SpeciesOfInterest == "TRUE")


# Get order of sp
Species_Order <- Freq_Taxa_min2_Samples %>%
  filter(Species %in% Freq_Taxa_min2_Samples_long_sub$Species) %>%
  arrange(Prevalence_URT)%>% 
  pull(Species) 

# Factor main df with this sp order
Freq_Taxa_min2_Samples_long_sub$Species <- factor(Freq_Taxa_min2_Samples_long_sub$Species, 
                                                     levels = Species_Order
                                                     ,ordered = T)
# plot
ggObj <- ggplot(Freq_Taxa_min2_Samples_long_sub, aes(x=SampleType, y=Species, fill = Prevalence_SampleType ))

savePlot <- ggObj + geom_tile() + theme_light(base_size = 14) + 
  labs(x="Sample Type", title = "Top 25 prevalent species in BAL", y= "Species") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        #axis.text.y=element_blank(), #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
  )

# save plot
ggsave(file="./Plots/SampleType_Prevalence_Top25_BALF_prevalent_taxa.png", plot=savePlot, width = 45, height = 20, units = "cm")

# Try boxplot
# plot
  ggObj <- ggplot(Freq_Taxa_min2_Samples_long, aes(x=SampleType, y =Prevalence_SampleType ))

savePlot <- ggObj + 
  geom_boxplot(alpha = 0.5) +
  geom_point() +
  labs(x="Sample Type", title = "Prevalence v/s Sampletype", y= "Top 25 Species prevalence in BAL") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=16,face="bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        #axis.text.y=element_blank(), #remove y axis labels
        #axis.ticks.y=element_blank()  #remove y axis ticks
  ) 


# save plot
ggsave(file="./Plots/Boxplot_SampleType_Prevalence_Top25_BALF_prevalent_taxa.png", plot=savePlot, width = 32, height = 15, units = "cm")










