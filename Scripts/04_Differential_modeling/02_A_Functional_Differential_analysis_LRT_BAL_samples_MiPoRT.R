# Script to detect KO associations with MaAsLin 3
# Mirrors the condensed LRT species workflow, but uses KO feature tables instead.
# Main idea:
#   1. Read KO feature table + metadata
#   2. Keep LRT samples only
#   3. Keep Healthy / Pneumonia / Covid-19 only
#   4. Remove the single-sample healthy project (PRJNA655567)
#   5. Create Dataset_Group
#   6. Prevalence-filter KOs
#   7. Run MaAsLin-3 separately for each disease dataset vs Healthy-Grp

# ------------------------------------------------------------------
# 0. Load libraries
# ------------------------------------------------------------------
for (lib in 
     c("maaslin3", "dplyr", "tidyr", "purrr", "ggplot2", "knitr")) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# ------------------------------------------------------------------
# 1. User-defined input paths and settings
# ------------------------------------------------------------------

# KO feature table
KO_TABLE_PATH <- "../01_KO_abundance_filtering/LRT_featureTable_KO_filt_AbundMeanAtLeastOneDisease_1eminus04.tsv"

# Metadata table
METADATA_PATH <- "../../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv"

# Name of the first column in the KO table that contains KO IDs
# Example possibilities: "KO", "KO_ext", "KO_ext_relab", "feature", etc.
FEATURE_ID_COL <- "Feature"

# Metadata sample ID column
SAMPLE_ID_COL <- "SampleID"

# Diversity column
DIVERSITY_COL <- "Richness"

# Output directory prefix
OUTPUT_PREFIX <- "LRT_KO"

# Number of CPU cores for MaAsLin3
N_CORES <- 6

# Reproducibility
set.seed(8284)

# ------------------------------------------------------------------
# 2. Read KO feature table
# ------------------------------------------------------------------
KO_profiles <- read.table(
    KO_TABLE_PATH,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

# Make feature IDs the row names
if (!(FEATURE_ID_COL %in% colnames(KO_profiles))) {
    stop(paste0("FEATURE_ID_COL '", FEATURE_ID_COL, "' not found in KO table."))
}

rownames(KO_profiles) <- KO_profiles[[FEATURE_ID_COL]]
KO_profiles[[FEATURE_ID_COL]] <- NULL

# Sanity checks
nSamples <- ncol(KO_profiles)
nFeatures <- nrow(KO_profiles)

print(paste0("Check if these are the number of samples you have: ", nSamples))
print(paste0("Check if these are the number of KO features you have: ", nFeatures))

# ------------------------------------------------------------------
# 3. Read metadata and keep LRT only
# ------------------------------------------------------------------
Metadata_Df <- read.table(
    METADATA_PATH,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

# Remove blank rows from metadata
Metadata_Df <- Metadata_Df[rowSums(is.na(Metadata_Df)) < ncol(Metadata_Df), ]

# Keep LRT only
Metadata_Df <- Metadata_Df %>%
    filter(RT_category == "LRT")

# Check that required columns exist
required_metadata_cols <- c(
    SAMPLE_ID_COL, "Disease", "BioProject", "After_QC_R1", "After_QC_R2", DIVERSITY_COL
)

missing_metadata_cols <- setdiff(required_metadata_cols, colnames(Metadata_Df))
if (length(missing_metadata_cols) > 0) {
    stop(
        paste(
            "Missing required metadata columns:",
            paste(missing_metadata_cols, collapse = ", ")
        )
    )
}

# Keep only diseases of interest
Disease_of_Interest <- c("Healthy", "Pneumonia", "Covid-19")

Metadata_Df_DOI <- Metadata_Df %>%
    filter(Disease %in% Disease_of_Interest) %>%
    filter(.data[[DIVERSITY_COL]] > 2)

# Add row names
rownames(Metadata_Df_DOI) <- Metadata_Df_DOI[[SAMPLE_ID_COL]]

# Quick checks
print(dim(Metadata_Df_DOI))
print(table(Metadata_Df_DOI$Disease, useNA = "ifany"))
print(table(Metadata_Df_DOI$BioProject, Metadata_Df_DOI$Disease, useNA = "ifany"))

# ------------------------------------------------------------------
# 4. Remove the single-sample project and define dataset groups
# ------------------------------------------------------------------
Metadata_Df_DOI_sub <- Metadata_Df_DOI %>%
    filter(BioProject != "PRJNA655567") %>%
    mutate(
        Total_reads = After_QC_R1 + After_QC_R2
    ) %>%
    select(all_of(c(SAMPLE_ID_COL, "BioProject", "Total_reads", "PatientID", "Disease"))) %>%
    mutate(
        Dataset_Group = case_when(
            BioProject == "PRJEB29011" ~ "Pneumonia-Dataset-1",
            BioProject == "PRJNA413615" ~ "Pneumonia-Dataset-2",
            BioProject == "PRJNA687506" ~ "Covid-19-Dataset",
            BioProject %in% c("PRJNA757846", "UniBergen") ~ "Healthy-Grp",
            TRUE ~ "Other"
        ),
        Dataset_Group = factor(Dataset_Group,
                               levels = c(
                                   "Healthy-Grp",
                                   "Covid-19-Dataset",
                                   "Pneumonia-Dataset-1",
                                   "Pneumonia-Dataset-2",
                                   "Other"
                                   )),
        Dataset_Group = droplevels.factor(Dataset_Group)
    )

print(table(Metadata_Df_DOI_sub$Dataset_Group, useNA = "ifany"))
"
        Healthy-Grp    Covid-19-Dataset Pneumonia-Dataset-1 
                 32                 105                 139 
Pneumonia-Dataset-2 
                 15
"
# The sample size looks consistent with the species level Maaslin3 models
# Hurray! 

print(dim(Metadata_Df_DOI_sub))

# ------------------------------------------------------------------
# 5. Match KO table to metadata samples
# ------------------------------------------------------------------
samples_present <- Metadata_Df_DOI_sub[[SAMPLE_ID_COL]] %in% colnames(KO_profiles)
if (!all(samples_present)) {
    missing_samples <- Metadata_Df_DOI_sub[[SAMPLE_ID_COL]][!samples_present]
    stop(
        paste(
            "These metadata samples are missing in the KO table:",
            paste(missing_samples, collapse = ", ")
        )
    )
}

KO_profiles_DOI <- KO_profiles %>%
    select(all_of(Metadata_Df_DOI_sub[[SAMPLE_ID_COL]]))

print(dim(KO_profiles_DOI))

# ------------------------------------------------------------------
# 6. Add factor levels for disease
# ------------------------------------------------------------------
#  sanity check
table(Metadata_Df_DOI_sub$Dataset_Group, 
      Metadata_Df_DOI_sub$Disease, 
      useNA = 'ifany')

# Add factors to Disease variable
Metadata_Df_DOI_sub$Disease <- factor(
    Metadata_Df_DOI_sub$Disease,
    levels = c("Healthy", "Covid-19", "Pneumonia")
)

Unique_Runs <- unique(Metadata_Df_DOI_sub$Dataset_Group)
print(Unique_Runs)

# We only want each disease dataset vs Healthy-Grp
Unique_Runs <- Unique_Runs[Unique_Runs != "Healthy-Grp" & Unique_Runs != "Other"]

# Drop unwanted levels
Unique_Runs <- droplevels.factor(Unique_Runs)
print(Unique_Runs)

# ------------------------------------------------------------------
# 7. Run MaAsLin3 separately per dataset group
# ------------------------------------------------------------------
for (each_Grp in Unique_Runs) {
    
    message("--------------------------------------------------")
    message("Running dataset group: ", each_Grp)
    
    # Keep one disease dataset + healthy group
    Metadata_sub <- Metadata_Df_DOI_sub %>%
        filter(Dataset_Group %in% c(each_Grp, "Healthy-Grp"))
    
    # Drop unused factor levels
    Metadata_sub$Disease <- droplevels(Metadata_sub$Disease)
    
    # Add row names
    rownames(Metadata_sub) <- Metadata_sub[[SAMPLE_ID_COL]]
    
    # Determine disease name for coef plotting
    Disease_Name <- levels(Metadata_sub$Disease)[2]
    
    # Filter KO table to:
    #   1. retained KOs
    #   2. samples in this metadata subset
    KO_profiles_sub <- KO_profiles %>%
        select(all_of(Metadata_sub[[SAMPLE_ID_COL]]))
    
    # Extra alignment check
    if (!identical(colnames(KO_profiles_sub), Metadata_sub[[SAMPLE_ID_COL]])) {
        KO_profiles_sub <- KO_profiles_sub[, Metadata_sub[[SAMPLE_ID_COL]], drop = FALSE]
    }
    
    print(dim(KO_profiles_sub))
    print(table(Metadata_sub$Dataset_Group, useNA = "ifany"))
    print(table(Metadata_sub$Disease, useNA = "ifany"))
    
    output_Name <- paste0(OUTPUT_PREFIX, "_", each_Grp, "_Model-2")
    
    RT_fit_out <- maaslin3(
        input_data = KO_profiles_sub,
        input_metadata = Metadata_sub,
        output = output_Name,
        formula = "~ Disease + Total_reads",
        normalization = "TSS",
        transform = "LOG",
        augment = TRUE,
        standardize = TRUE,
        max_significance = 0.1,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = TRUE,
        max_pngs = 500,
        cores = N_CORES,
        coef_plot_vars = c(paste0("Disease ", Disease_Name)),
        heatmap_vars = c("Total_reads")
    )
}

# ------------------------------------------------------------------
# 8. Run MaAsLin3 for all Data
# ------------------------------------------------------------------

output_Name <- paste0(OUTPUT_PREFIX, "_", "Model-1")

# sanity check
table(
    Metadata_Df_DOI_sub$BioProject,
    Metadata_Df_DOI_sub$Disease,
    useNA = 'ifany'
      )

levels(Metadata_Df_DOI_sub$Disease)

# Fit Global model
RT_fit_out <- maaslin3(
    input_data = KO_profiles_DOI,
    input_metadata = Metadata_Df_DOI_sub,
    output = output_Name,
    formula = "~ Disease + Total_reads",
    normalization = "TSS",
    transform = "LOG",
    augment = TRUE,
    standardize = TRUE,
    max_significance = 0.1,
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = TRUE,
    max_pngs = 500,
    cores = N_CORES,
    coef_plot_vars = c("Disease Pneumonia",
                       "Disease Covid-19"),
    heatmap_vars = c("Total_reads")
)

# ------------------------------------------------------------------
# 9. End
# ------------------------------------------------------------------
message("Done running KO-based LRT MaAsLin3 models.")
