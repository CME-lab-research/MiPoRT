# Script to filter merged feature table with KO profiles
# Specifically, this script would 
# 1. remove UNGROUPED
# 2. remove zero-sum samples
# 3. remove all-zero features
# 4. remove samples which failed QC checks (using Metadata table)
# 5. remove features present in <10 samples


library(tidyverse)

# ------------- # ------------- # ------------- # ------------- #
# ------------- # ------------- # ------------- # ------------- #
#                            Part A:                            #
# ------------- # ------------- # ------------- # ------------- #
# ------------- # ------------- # ------------- # ------------- #

# Part A: Filtering standard KO mapping file

# Step-1 : Read merged HUMANN pathway abundance table
# NOTE: Rm the '#' symbol from the Merged Humann file
KO_raw <- read.table("./03_A_Humann_merged_table/Humann4_merged_KO.tsv",
                     sep = '\t',
                     header = T)

# Step 1: drop the 'UNGROUPED' feature
    KO_raw <- KO_raw[KO_raw$Adjusted_CPMs != 'UNGROUPED',]
    row.names(KO_raw) <- KO_raw$Adjusted_CPMs
        
# logging and sanity check
writeLines(paste("Sample count:", dim(KO_raw)[2],"\nFeature count:", dim(KO_raw)[1]))
    # Sample count: 3504 
    # Feature count: 10663

# Step 2: remove zero-sum samples
    # Check number of samples with total abundance = 0
    table(colSums(KO_raw[,-1]) == 0, useNA = 'ifany')    
    "
       FALSE  TRUE 
        3473    30
    "
    # visual
    savePlot <- hist(colSums(KO_raw[,-1]), 
         breaks = 50, main = "Sample total abundances before filtering",
         xlab = "Total abundance", col = "skyblue")
    abline(v = 0.001, col = "red", lty = 2)

    # get names of samples to remove
    Samples_to_rm <- names(which((colSums(KO_raw[,-1]) == 0) == TRUE))

    # Remove 0 abundance i.e, empty samples
    KO_filt <- KO_raw %>%
        select(!all_of(c(Samples_to_rm)))
    
    # logging and sanity check
    writeLines(paste("Sample count:", dim(KO_filt)[2], 
                     "\nFeature count:", dim(KO_filt)[1]))
    # Sample count: 3474 
    # Feature count: 10663
    
    table(colSums(KO_filt[,-1]) == 0, useNA = 'ifany')  
    "
       FALSE
        3473
    "

# Step 3: remove all-zero features
    # Check how many are all zero
    table(rowSums(KO_raw[,-1]) == 0, useNA = 'ifany')    
    "
    FALSE
    10663
    "
# Step 4: Remove samples which failed QC checks (using Metadata table)
    
    # read Metadata table
    Metadata_df <- read.table(
        "./S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", header = T, sep = "\t")
    
    # check how many samples are in both tables
    RetainSample <- (colnames(KO_filt)[-1] %in% Metadata_df$SampleID)
    
    # logging
    table(RetainSample)
    # RetainSample
    # FALSE  TRUE 
    #   338  3135
    
    # filter
    KO_filt <- KO_filt %>%
        select(all_of(c("Adjusted_CPMs", Metadata_df$SampleID)))
    
    # logging and sanity check
    writeLines(paste("Sample count:", dim(KO_filt)[2]-1, 
                     "\nFeature count:", dim(KO_filt)[1]))
    # Sample count: 3135 
    # Feature count: 10663
    
    # visual
    hist(colSums(KO_filt[,-1]), 
         breaks = 50,
         main = "Sample total abundances after filtering",
         xlab = "Total abundance", col = "skyblue")
    abline(v = 0.001, col = "red", lty = 2)
    
# 5. remove features present in <10 samples
    row.names(KO_filt) <- KO_filt$Adjusted_CPMs
    KO_filt <- KO_filt[,-1]
    
    # check total KO prevalence
    KO_terms_Total_prevalence <- rowSums(KO_filt > 0)
    KO_terms_Total_abundance <- rowSums(KO_filt)
    
    # Find number of species with at least x total abundance
    hist(KO_terms_Total_abundance/3135,
         breaks = 50,
         main = "KO total abundances after filtering",
         xlab = "Total abundance", col = "skyblue")    
    
    # Features passing 0.1% sample count i.e, (N=3 samples)
    table(KO_terms_Total_prevalence/3135 >= 0.001, useNA = 'ifany')
    # 9181 features

    # Features passing 0.5% sample count i.e, (N=15 samples)
    table(KO_terms_Total_prevalence/3135 >= 0.005, useNA = 'ifany')
    # 7787 features

    # Features passing 1% sample count i.e, (N=31 samples)
    table(KO_terms_Total_prevalence/3135 >= 0.01, useNA = 'ifany')
    # 7056 features
    
    # Abundance filter checks
    
    # Features passing 0.001 mean abundance
    table(KO_terms_Total_abundance/3135 >= 0.001, useNA = 'ifany')
    # 9698 features
    
    # Features passing 0.005 mean abundance
    table(KO_terms_Total_abundance/3135 >= 0.005, useNA = 'ifany')
    # 9185 features

    # Features passing 0.05 mean abundance    
    table(KO_terms_Total_abundance/3135 >= 0.05, useNA = 'ifany')
    # 7644 features    
    
    # Mean abundance can be misleading. Filter based on prevalence.
    # Filter that combines prevalence and abundance cutoffs
    
    # Filter out species with <1 CPM abundance per 10 samples
    min_samples = 10
    KO_prevalence_filtered <- KO_filt[
        rowSums(KO_filt >= 1) >= min_samples, 
    ]
    
    # keeps species present in at least 10 samples (min_samples) with a  abundance >= 1 (Adjusted CPMs)
    
    # logging and sanity check
    writeLines(paste("Sample count:", dim(KO_prevalence_filtered)[2], 
                     "\nFeature count:", dim(KO_prevalence_filtered)[1]))
    # Sample count: 3135 
    # Feature count: 5601


# Remove 0 variance features
    # remove zero-variance features (important)
    KO_var <- apply(KO_prevalence_filtered, 1, var)
    KO_prevalence_filtered <- KO_prevalence_filtered[KO_var > 0, ]

    # logging and sanity check
    writeLines(paste("Sample count:", dim(KO_prevalence_filtered)[2], 
                     "\nFeature count:", dim(KO_prevalence_filtered)[1]))
    # Sample count: 3135 
    # Feature count: 5601
    
    # Extra visual check
    # Compute sample sums
    df_plot <- data.frame(
        Sample_Total = c(
            colSums(KO_raw[, -1]),
            colSums(KO_prevalence_filtered[,])),
        Dataset = rep(c("Before filtering", "After filtering"),
                      times = c(
                          ncol(KO_raw[, -1]),
                          ncol(KO_prevalence_filtered[,])))
    )
    
    # add factors
    # sanity check
    levels(df_plot$Dataset)
    
    # mutate
    df_plot <- df_plot %>%
        mutate(Dataset = factor(Dataset,
                                levels = c("Before filtering", "After filtering")))
    
    # sanity check
    levels(df_plot$Dataset)
    
    # Plot
    savePlot <- ggplot(df_plot, aes(x = Sample_Total)) +
        geom_histogram(bins = 50, fill = "skyblue", color = "black") +
        facet_wrap(~ Dataset, scales = "fixed") +
        labs(
            title = "Sample total KO abundances",
            x = "Total abundance",
            y = "Count"
        ) +
        theme_bw() + 
        theme(axis.text=element_text(size=20), 
              axis.title=element_text(size=20,face="bold"),
              title = element_text(size=20,face="bold"),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 18, face = "bold"),
              panel.spacing=unit(2,"cm"),
              plot.margin = margin(t = 20, 
                                   r = 40, 
                                   b = 20, 
                                   l = 40, 
                                   unit = "pt"),
              strip.text.x = element_text(size = 18, face = "bold"),
              strip.background = element_blank()
              ) +
        scale_x_log10() 
    
    savePlot
    
    # Save the plot
    ggsave("./Plots/Histogram_Total_KO_abundances_rawCounts.png", 
           plot=savePlot, 
           height = 16, 
           width = 36, 
           units = "cm", 
           limitsize = FALSE, 
           dpi = 150)
    
    # Save the plot log 10
    ggsave("./Plots/Histogram_Total_KO_abundances_Log10_Counts.png", 
           plot=savePlot, 
           height = 16, 
           width = 36, 
           units = "cm", 
           limitsize = FALSE, 
           dpi = 150)

    # logging and sanity check
    writeLines(paste("Sample count:", dim(KO_prevalence_filtered)[2], 
                     "\nFeature count:", dim(KO_prevalence_filtered)[1]))
    # Sample count: 3135 
    # Feature count: 5601
    
    
    # save the table and proceed with normalization and filtering
write.table(KO_prevalence_filtered,
            "./01_Humann4_merged_KO_filtered.tsv",
            row.names = T,
            sep = '\t',
            quote = FALSE)

# ------------- # ------------- # ------------- # ------------- #
# ------------- # ------------- # ------------- # ------------- #
#                            Part B:                            #
# ------------- # ------------- # ------------- # ------------- #
# ------------- # ------------- # ------------- # ------------- #

# Part B: Filtering extended KO mapping based feature table

# Open the file and remove comment char

# Step-1 : Read merged HUMANN pathway abundance table
# NOTE: Rm the '#' symbol from the Merged Humann file
KO_raw <- read.table("./03_A_Humann_merged_table/Humann4_merged_KO_expanded.tsv",
                     sep = '\t',
                     header = T)

# Step 1: drop the 'UNGROUPED' feature
KO_raw <- KO_raw[KO_raw$Adjusted_CPMs != 'UNGROUPED',]
row.names(KO_raw) <- KO_raw$Adjusted_CPMs
KO_raw <- KO_raw[,-1]

# logging and sanity check
writeLines(paste("Sample count:", dim(KO_raw)[2],"\nFeature count:", dim(KO_raw)[1]))
# Sample count: 3503 
# Feature count: 11949

# Step 2: remove zero-sum samples
# Check number of samples with total abundance = 0
table(colSums(KO_raw) == 0, useNA = 'ifany')    
"
       FALSE  TRUE 
        3473    23
    "
# visual
savePlot <- hist(colSums(KO_raw), 
                 breaks = 50, main = "Sample total abundances before filtering",
                 xlab = "Total abundance", col = "skyblue")
abline(v = 0.001, col = "red", lty = 2)

# get names of samples to remove
Samples_to_rm <- names(which((colSums(KO_raw) == 0) == TRUE))

# Remove 0 abundance i.e, empty samples
KO_filt <- KO_raw %>%
    select(!all_of(c(Samples_to_rm)))

# logging and sanity check
writeLines(paste("Sample count:", dim(KO_filt)[2], 
                 "\nFeature count:", dim(KO_filt)[1]))
# Sample count: 3480 
# Feature count: 11949

table(colSums(KO_filt) == 0, useNA = 'ifany')  
"
       FALSE
        3480
    "

# Step 3: remove all-zero features
# Check how many are all zero
table(rowSums(KO_raw) == 0, useNA = 'ifany')    
"
    FALSE
    11949
"

# Step 4: Remove samples which failed QC checks (using Metadata table)

# read Metadata table
Metadata_df <- read.table(
    "./S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv", header = T, sep = "\t")

# check how many samples are in both tables
RetainSample <- (colnames(KO_filt) %in% Metadata_df$SampleID)

# logging
table(RetainSample)
# RetainSample
# FALSE  TRUE 
#   345  3135

# filter
KO_filt <- KO_filt %>%
    select(all_of(c(Metadata_df$SampleID)))

# logging and sanity check
writeLines(paste("Sample count:", dim(KO_filt)[2]-1, 
                 "\nFeature count:", dim(KO_filt)[1]))

# Sample count: 3134 
# Feature count: 11949

# visual
hist(colSums(KO_filt), 
     breaks = 50,
     main = "Sample total abundances after filtering",
     xlab = "Total abundance", col = "skyblue")
abline(v = 0.001, col = "red", lty = 2)

# 5. remove features present in <10 samples
    # check total KO prevalence
    KO_terms_Total_prevalence <- rowSums(KO_filt > 0)
    KO_terms_Total_abundance <- rowSums(KO_filt)

# Find number of species with at least x total abundance
hist(KO_terms_Total_abundance/3135,
     breaks = 50,
     main = "KO total abundances after filtering",
     xlab = "Total abundance", col = "skyblue")    

    # Features passing 0.1% sample count i.e, (N=3 samples)
    table(KO_terms_Total_prevalence/3135 >= 0.001, useNA = 'ifany')
    # 10375 features
    
    # Features passing 0.5% sample count i.e, (N=15 samples)
    table(KO_terms_Total_prevalence/3135 >= 0.005, useNA = 'ifany')
    # 8988 features
    
    # Features passing 1% sample count i.e, (N=31 samples)
    table(KO_terms_Total_prevalence/3135 >= 0.01, useNA = 'ifany')
    # 8256 features
    
    # Abundance filter checks
    
    # Features passing 0.001 mean abundance
    table(KO_terms_Total_abundance/3135 >= 0.001, useNA = 'ifany')
    # 10882 features
    
    # Features passing 0.005 mean abundance
    table(KO_terms_Total_abundance/3135 >= 0.005, useNA = 'ifany')
    # 10348 features    

    # Features passing 0.05 mean abundance    
    table(KO_terms_Total_abundance/3135 >= 0.05, useNA = 'ifany')
    # 8924 features    

# Mean abundance can be misleading. Filter based on prevalence.
# Filter that combines prevalence and abundance cutoffs

    # Filter out species with <1 CPM abundance per 10 samples
    min_samples = 10
    KO_prevalence_filtered <- KO_filt[
        rowSums(KO_filt >= 1) >= min_samples, 
    ]

    # keeps species present in at least 10 samples (min_samples) with a  abundance >= 1 (Adjusted CPMs)

# logging and sanity check
writeLines(paste("Sample count:", dim(KO_prevalence_filtered)[2], 
                 "\nFeature count:", dim(KO_prevalence_filtered)[1]))
# Sample count: 3135 
# Feature count: 9048

# Remove 0 variance features
# remove zero-variance features (important)
KO_var <- apply(KO_prevalence_filtered, 1, var)
KO_prevalence_filtered <- KO_prevalence_filtered[KO_var > 0, ]

# logging and sanity check
writeLines(paste("Sample count:", dim(KO_prevalence_filtered)[2], 
                 "\nFeature count:", dim(KO_prevalence_filtered)[1]))
# Sample count: 3135 
# Feature count: 9048

# Extra visual check
# Compute sample sums
df_plot <- data.frame(
    Sample_Total = c(
        colSums(KO_raw),
        colSums(KO_prevalence_filtered)),
    Dataset = rep(c("Before filtering", "After filtering"),
                  times = c(
                      ncol(KO_raw),
                      ncol(KO_prevalence_filtered)))
)

# add factors
# sanity check
levels(df_plot$Dataset)

# mutate
df_plot <- df_plot %>%
    mutate(Dataset = factor(Dataset,
                            levels = c("Before filtering", "After filtering")))

# sanity check
levels(df_plot$Dataset)

# Plot
savePlot <- ggplot(df_plot, aes(x = Sample_Total)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black") +
    facet_wrap(~ Dataset, scales = "fixed") +
    labs(
        title = "Sample total KO (extended) abundances",
        x = "Total abundance",
        y = "Count"
    ) +
    theme_bw() + 
    theme(axis.text=element_text(size=20), 
          axis.title=element_text(size=20,face="bold"),
          title = element_text(size=20,face="bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold"),
          panel.spacing=unit(2,"cm"),
          plot.margin = margin(t = 20, 
                               r = 40, 
                               b = 20, 
                               l = 40, 
                               unit = "pt"),
          strip.text.x = element_text(size = 18, face = "bold"),
          strip.background = element_blank()
    ) +
    scale_x_log10() 

savePlot

# Save the plot
ggsave("./Plots/Histogram_Total_KO_extended_abundances_rawCounts.png", 
       plot=savePlot, 
       height = 16, 
       width = 36, 
       units = "cm", 
       limitsize = FALSE, 
       dpi = 150)

# Save the plot log 10
ggsave("./Plots/Histogram_Total_KO_extended_abundances_Log10_Counts.png", 
       plot=savePlot, 
       height = 16, 
       width = 36, 
       units = "cm", 
       limitsize = FALSE, 
       dpi = 150)

# logging and sanity check
writeLines(paste("Sample count:", dim(KO_prevalence_filtered)[2], 
                 "\nFeature count:", dim(KO_prevalence_filtered)[1]))
# Sample count: 3135 
# Feature count: 9048


# save the table and proceed with normalization and filtering
write.table(KO_prevalence_filtered,
            "./01_Humann4_merged_KO_extended_filtered.tsv",
            row.names = T,
            sep = '\t',
            quote = FALSE)


