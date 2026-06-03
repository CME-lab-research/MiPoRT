# ============================================================
# 0. Load libraries
# ============================================================

# ---- 0.1 Load required packages ----

# Load core tidyverse-style packages
for (lib in c(
  "dplyr",
  "scales",
  "tidyr",
  "purrr",
  "ggplot2",
  "tibble",
  "pheatmap",
  "RColorBrewer",
  "ghibli",
  "forcats"
)) {
  suppressPackageStartupMessages(
    require(lib, character.only = TRUE)
  )
}

# ---- 0.2 Load custom metadata factor function ----

# This function adds ordered factors for:
# - SampleType
# - RT_category
# - AgeGroup
# etc.

source("../Scripts/add_metadata_factors.R")


# ============================================================
# 1. User settings
# ============================================================

# ---- 1.1 Input files ----

# KO abundance table from HUMAnN
KO_TABLE_PATH <- "../03_C_Humann_BatchCorr/MiPORT_BC-RT_withAN_featureTable_KO_ext.tsv"

# Metadata file
METADATA_PATH <- "../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv"

# KO to KEGG BRITE mapping table
KO_BRITE_MAPPING_PATH <- "C:/Users/o_shinde/Desktop/KO_BRITE_map/ko00001_KO_mapping.tsv"

# ---- 1.2 Output directory ----

OUTDIR <- "Plots"
dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)

# ---- 1.3 Define feature ID column ----

FEATURE_ID_COL <- "KO_ext_relab"

# ---- 1.4 Define metadata sample ID column ----

SAMPLE_ID_COL <- "SampleID"

# ---- 1.5 SampleTypes to include ----

SampleTypes_Of_Interest <- c(
  "Nasal_Swab",
  "Buccal_mucosa",
  "Oral_swab",
  "Saliva",
  "Tongue_dorsum",
  "Sputum",
  "BAL"
)


# ============================================================
# 2. Read KO feature table
# ============================================================

# ---- 2.1 Read HUMAnN KO abundance table ----

KO_profiles <- read.table(
  KO_TABLE_PATH,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# ---- 2.2 Sanity check feature column ----

if (!(FEATURE_ID_COL %in% colnames(KO_profiles))) {
  
  stop(
    paste0(
      "Feature ID column not found: ",
      FEATURE_ID_COL
    )
  )
}

# ---- 2.3 Set KO IDs as rownames ----

rownames(KO_profiles) <- KO_profiles[[FEATURE_ID_COL]]

# Remove feature column after converting to rownames
KO_profiles[[FEATURE_ID_COL]] <- NULL

# ---- 2.4 Logging ----

cat("KO table loaded.\n")
cat("Number of KO features:", nrow(KO_profiles), "\n")
cat("Number of samples:", ncol(KO_profiles), "\n\n")


# ============================================================
# 3. Read metadata
# ============================================================

# ---- 3.1 Read metadata ----

Metadata_Df <- read.table(
  METADATA_PATH,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# ---- 3.2 Remove empty rows ----

Metadata_Df <- Metadata_Df[
  rowSums(is.na(Metadata_Df)) < ncol(Metadata_Df),
]

# ---- 3.3 Add metadata factors ----

Metadata_Df <- add_metadata_factors(Metadata_Df)

# ---- 3.4 Define required metadata columns ----

required_metadata_cols <- c(
  "BioProject",
  "Citation",
  SAMPLE_ID_COL,
  "SampleTypev2",
  "RT_category",
  "Disease"
)

# ---- 3.5 Sanity check required columns ----

missing_cols <- setdiff(
  required_metadata_cols,
  colnames(Metadata_Df)
)

if (length(missing_cols) > 0) {
  
  stop(
    paste(
      "Missing required metadata columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
  
} else {
  
  Metadata_Df <- Metadata_Df %>%
    select(all_of(required_metadata_cols))
}

# ---- 3.6 Filter metadata ----

# Keep only:
# - Healthy samples
# - select respiratory sample types
Metadata_Df_Healthy <- Metadata_Df %>%
  filter(
    Disease == "Healthy",
    SampleTypev2 %in% SampleTypes_Of_Interest
  )

  # sanity check
  Metadata_Df_Healthy %>% 
    count(SampleTypev2)
  '
     SampleTypev2   n
  1           BAL  50
  2 Buccal_mucosa 527
  3    Nasal_Swab  18
  4     Oral_swab 171
  5        Saliva 210
  6        Sputum  15
  7 Tongue_dorsum 418
  '
  
# ---- 3.7 Logging ----

cat("Metadata loaded.\n")
cat(
  "Healthy samples retained:",
  nrow(Metadata_Df_Healthy),
  "\n\n"
)


# ============================================================
# 4. Read KO-BRITE mapping
# ============================================================

# ---- 4.1 Read KO mapping file ----

KO_BRITE_map <- read.table(
  KO_BRITE_MAPPING_PATH,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# ---- 4.2 Expected columns ----

# Expected example:
# KO
# KO_Name
# KEGG_category
# KEGG_subcategory

required_mapping_cols <- c(
  "KO",
  "KEGG_category_name",
  "KEGG_subcategory_name"
)

missing_mapping_cols <- setdiff(
  required_mapping_cols,
  colnames(KO_BRITE_map)
)

if (length(missing_mapping_cols) > 0) {
  
  stop(
    paste(
      "Missing mapping columns:",
      paste(missing_mapping_cols, collapse = ", ")
    )
  )
}

# ---- 4.2 Clean KO-BRITE mapping table ----

# List of KO categories to keep
KO_BRITE_A_to_Keep <- c(
  "Metabolism",
  "Organismal Systems",
  "Cellular Processes",
  "Genetic Information Processing",
  "Environmental Information Processing",
  "Human Diseases"
)

# Filter and clean
KO_BRITE_map_clean <- 
  KO_BRITE_map %>%
  select(
    KO,
    gene_symbol,
    description,
    KEGG_category_name,
    KEGG_subcategory_name,
    KEGG_pathway_name
  ) %>%
  
  # Remove incomplete annotations
  filter(
    !is.na(KO),
    !is.na(KEGG_subcategory_name),
    !is.na(KEGG_pathway_name)
  ) %>%
  
  # filter
  filter(KEGG_category_name %in% KO_BRITE_A_to_Keep) %>%

  # Remove exact duplicate KO-pathway-category rows
  distinct()

# ---- 4.3 Logging ----

cat("KO-BRITE mapping loaded.\n")
cat(
  "Mapped KOs:",
  nrow(KO_BRITE_map_clean),
  "\n\n"
)

  # Mapped KOs: 37980

# ============================================================
# 5. Subset KO table to filtered metadata samples
# ============================================================

# ---- 5.1 Identify overlapping samples ----

Samples_Keep <- intersect(
  colnames(KO_profiles),
  Metadata_Df_Healthy[[SAMPLE_ID_COL]]
)

# ---- 5.2 Subset KO table ----

# Subset for sampleID's of interest
KO_profiles_sub <- KO_profiles[
  ,
  Samples_Keep,
  drop = FALSE
]

# Subset for KO's of interest
# Select top KOs (Prevalent >20% in any sampleType)
# Read prevalence calculation file to get a list of KOs to keep
  Ko_prevalence_list <-
  read.table("./SampleType_KO_prevalence_stats.txt",
             header = T,
             sep = '\t')

  # Get unique KOs
  cat(
    "Number of KOs to retain:",
    length(unique(Ko_prevalence_list$Feature)),
    "\n\n"
  )
  # Before filtering for classes of interest
  # Number of KOs to retain: 4855

  KOs_to_retain <- unique(Ko_prevalence_list$Feature)
  
  # remove KOs not present in the clean KO_BRITE map
  KOs_to_retain <- c(intersect(KOs_to_retain,
                         unique(KO_BRITE_map_clean$KO)))
  
  # Get unique KOs
  cat(
    "Number of KOs retained after cleaning:",
    length(KOs_to_retain),
    "\n\n"
  )
  
  # Number of KOs retained after cleaning: 2763
  
  # subset
  KO_profiles_sub <- KO_profiles_sub %>%
    filter(row.names(KO_profiles_sub) %in% KOs_to_retain)
  
# ---- 5.3 Logging ----
  
  cat(
    "Samples retained in KO table:",
    ncol(KO_profiles_sub),
    "\nKOs retained in KO table:",
    nrow(KO_profiles_sub),
    "\n\n"
  )
  
  # Samples retained in KO table: 1409 
  # KOs retained in KO table: 2763
  
# ============================================================
# 6. Convert KO table to long format
# ============================================================

# ---- 6.1 Convert rownames to column ----
KO_profiles_long <- KO_profiles_sub %>%
  rownames_to_column("KO")

# ---- 6.2 Pivot longer ----
KO_profiles_long <- KO_profiles_long %>%
  pivot_longer(
    cols = -KO,
    names_to = "SampleID",
    values_to = "KO_abundance"
  )
  # should ideally give you 2763*1409 = 3893067 rows
  
# ---- 6.3 Join metadata ----

KO_profiles_long <- KO_profiles_long %>%
  left_join(
    Metadata_Df_Healthy,
    by = c("SampleID" = SAMPLE_ID_COL)
  )


# ---- 6.4 Logging ----

cat(
  "Rows in long-format table:",
  nrow(KO_profiles_long),
  "\n\n"
)
  
  # Rows in long-format table: 3893067


# ============================================================
# 7. Calculate mean KO abundance per SampleType
# ============================================================

# ---- 7.1 Calculate mean abundance ----

KO_abundance_mean_perSampleType <- 
    KO_profiles_long %>%
    group_by(
      SampleTypev2,
      RT_category,
      KO
    ) %>%
    summarise(
      KO_mean_abundance = mean(
        KO_abundance,
        na.rm = TRUE
      ),
      nSampleTypes = n_distinct(SampleID),
      .groups = "drop"
    )

#---- 7.2 Join KO-BRITE mapping ----
  
  # This is intentionally a many-to-many join.
  # A KO can belong to multiple KEGG pathways.
  # We keep all valid KO-pathway memberships here.
  KO_abundance_with_BRITE <- 
    KO_abundance_mean_perSampleType %>%
    inner_join(
      KO_BRITE_map_clean,
      by = "KO",
      relationship = "many-to-many"
    )

# [Intermediate summary] Calculate gene summary per SampleType
# ---- 7.3 Calculate mean KO abundance per SampleType
  # For each SampleType:
  # Gene_mean = mean abundance across SampleType replicates
  KEGG_gene_summary <- KO_abundance_with_BRITE %>%
    group_by(
      SampleTypev2,
      nSampleTypes,
      RT_category,
      KO,
      description,
      KEGG_category_name,
      KEGG_subcategory_name,
      KEGG_pathway_name
    ) %>%
    summarise(
      gene_mean = mean(KO_mean_abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # reorder columns
    select(all_of(
      c(
        "RT_category",
        "SampleTypev2",
        "KEGG_category_name",
        "KEGG_subcategory_name",
        "KEGG_pathway_name",
        "KO",
        "description",
        "gene_mean",
        "nSampleTypes"
      )
    ))

  # ---- 7.3.1 Save SampleType-level KO summary table ----
  
  # Save this intermediate table for biological interpretation,
  # figure checking, and manuscript writing.
  write_tsv(
    KEGG_gene_summary,
    "KEGG_gene_summary_SampleType_level.tsv"
  )
  
# ============================================================
# 8. Calculate KEGG pathway abundance per SampleType
# ============================================================

# ---- 8.1 Calculate mean KO abundance within each KEGG pathway
# For each SampleType and KEGG pathway:
# pathway_mean = mean abundance of all member KOs
#
# This prevents larger pathways from automatically having larger values
# only because they contain more KOs.
  
  KEGG_pathway_summary <- KO_abundance_with_BRITE %>%
    group_by(
      SampleTypev2,
      RT_category,
      KEGG_category_name,
      KEGG_subcategory_name,
      KEGG_pathway_name,
      nSampleTypes
    ) %>%
    summarise(
      pathway_mean = mean(KO_mean_abundance, na.rm = TRUE),
      pathway_total = sum(KO_mean_abundance, na.rm = TRUE),
      n_KOs_in_pathway = n_distinct(KO),
      .groups = "drop"
    ) %>%
    # reorder columns
    select(all_of(
      c(
        "RT_category",
        "SampleTypev2",
        "nSampleTypes",
        "KEGG_category_name",
        "KEGG_subcategory_name",
        "KEGG_pathway_name",
        "pathway_mean",
        "pathway_total",
        "n_KOs_in_pathway"
        
      )
    ))

  # 8.1.1 save pathway-level summary table per sampleType ----
  
  # Save this intermediate table for biological interpretation,
  # figure checking, and manuscript writing.
  write_tsv(
    KEGG_pathway_summary,
    "KEGG_pathway_summary_SampleType_level.tsv"
  )
  
  
# ============================================================
# 9. Calculate BRITE Level 2 abundance per SampleType
# ============================================================
  
  # ---- 9.1 Average pathway means within each BRITE Level 2 category ----
  
  # For each SampleType and BRITE Level 2 category:
  # level2_category_abundance = mean(pathway_mean)
  #
  # This gives each KEGG pathway equal weight within a BRITE Level 2 category.
  BRITE_level2_summary <- KEGG_pathway_summary %>%
    group_by(
      SampleTypev2,
      nSampleTypes,
      RT_category,
      KEGG_category_name,
      KEGG_subcategory_name
    ) %>%
    summarise(
      level2_category_abundance = mean(pathway_mean, 
                                       na.rm = TRUE),
      level2_category_total_pathway_mean = sum(pathway_mean, 
                                               na.rm = TRUE),
      n_pathways_in_level2 = n_distinct(KEGG_pathway_name),
      n_KOs_total_in_level2 = sum(n_KOs_in_pathway, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # reorder columns
    select(all_of(
      c(
        "RT_category",
        "SampleTypev2",
        "nSampleTypes",
        "KEGG_category_name",
        "KEGG_subcategory_name",
        "level2_category_abundance", # mean abundance(pathway_mean)
        "level2_category_total_pathway_mean", # sum(pathway_mean)
        "n_pathways_in_level2",
        "n_KOs_total_in_level2"
        
      )
    ))
  
  # Save this intermediate table for biological interpretation,
  # figure checking, and manuscript writing.
  write_tsv(
    BRITE_level2_summary,
    "BRITE_level2_summary_SampleType_level.tsv"
  )
  
  
# ============================================================
# 10. Prepare plotting dataframe
# ============================================================

  # Make a plot_df copy
  BRITE_level2_plot_df <- BRITE_level2_summary 
  
  # ---- 10.2 Add SampleType factor ----
  # sanity check
  levels(BRITE_level2_plot_df$SampleTypev2)
  
  # init factors and labels
  SamplingSite_Factor <- c(
    "Nasal_Swab", 
    "Buccal_mucosa", 
    "Oral_swab", 
    "Saliva", 
    "Tongue_dorsum", 
    "Sputum", 
    "BAL")
  
  # Update labels
  SamplingSite_labels <- c(
    "Nasal swab", 
    "Buccal mucosa", 
    "Oral swab", 
    "Saliva", 
    "Tongue dorsum", 
    "Sputum", 
    "BAL")
  
  # add factors
  BRITE_level2_plot_df$SampleTypev2 <- 
    factor(BRITE_level2_plot_df$SampleTypev2, 
           levels = SamplingSite_Factor, 
           labels = SamplingSite_labels)  
  
  # sanity check
  levels(BRITE_level2_plot_df$SampleTypev2)
  
  # Add RT category factors
  # sanity check
  levels(BRITE_level2_plot_df$RT_category)
  
  # Add KEGG category levels
  BRITE_level2_plot_df$KEGG_category_name <- 
    factor(BRITE_level2_plot_df$KEGG_category_name,
           levels = KO_BRITE_A_to_Keep
                                   )
  BRITE_level2_plot_df %>%
    count(KEGG_category_name)
  
  levels(BRITE_level2_plot_df$KEGG_category_name)
  
  # update the plotting df
  BRITE_level2_plot_df <- BRITE_level2_plot_df %>%
    mutate(
      # This column controls the tile fill intensity.
      # It represents the total abundance signal for the BRITE level 2 category.
      level2_total_abundance_plot = level2_category_total_pathway_mean,
      
      # This column controls bubble size.
      # It represents the mean abundance within the BRITE level 2 category.
      level2_mean_abundance_plot = level2_category_abundance
    )
  
  
# ---- 10.3 Reorder KEGG subcategories within each KEGG category ----
  # Order BRITE categories by mean abundance across all SampleTypes.
  # Calculate one total abundance score per KEGG subcategory.
  # This score is summed across all sample types.
  # It will be used only for ordering the y-axis.
  KEGG_subcategory_order_df <- BRITE_level2_plot_df %>%
    group_by(
      KEGG_category_name,
      KEGG_subcategory_name
    ) %>%
    summarise(
      subcategory_total_abundance = sum(
        level2_total_abundance_plot,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    
    # Sort within each KEGG category.
    # Higher-abundance subcategories will appear closer to the top after factor reversal.
    arrange(
      KEGG_category_name,
      subcategory_total_abundance
    ) %>%
    
    # Create a unique plotting label.
    # This avoids problems if the same subcategory name appears in more than one KEGG category.
    mutate(
      KEGG_subcategory_plot_id = paste(
        KEGG_category_name,
        KEGG_subcategory_name,
        sep = "___"
      )
    )
  

  # Join the ordering information back to the main plotting dataframe.
  BRITE_level2_plot_df_ordered <- BRITE_level2_plot_df %>%
    left_join(
      KEGG_subcategory_order_df,
      by = c(
        "KEGG_category_name",
        "KEGG_subcategory_name"
      )
    ) %>%
    
    # Convert the unique plotting ID to a factor.
    # The factor levels determine the y-axis order.
    mutate(
      KEGG_subcategory_plot_id = factor(
        KEGG_subcategory_plot_id,
        levels = KEGG_subcategory_order_df$KEGG_subcategory_plot_id
      )
    )
  
  
# ============================================================
# 11. Plot BRITE heatmap-bubble plot
# ============================================================

# ---- 11.1 Generate ggplot ----
# Initialize the plot with sample type on x-axis
# and KEGG BRITE level 2 category on y-axis.

BRITE_bubble_plot <- ggplot(
    BRITE_level2_plot_df_ordered,
    aes(
      x = SampleTypev2,
      y = KEGG_subcategory_plot_id
    )
    )
  
  # build aesthetic details
  BRITE_bubble_plot <- 
    BRITE_bubble_plot +
  
    # Draw heatmap-style tiles.
    # Tile fill shows the total abundance of pathways/KOs
    # summarized into each BRITE level 2 category.
    geom_tile(
      aes(fill = level2_total_abundance_plot),
      color = "white",
      linewidth = 0.4
    ) +
    
    # Use a white-to-magenta gradient for summed mean relative abundance.
    # The values are displayed as percentages because the original data
    # came from relative abundance profiles.
    scale_fill_gradient(
      low = "#FFF5FB",
      high = "#8B008B",
      trans = pseudo_log_trans(base = 10),
      limits = c(
        0,
        quantile(
          BRITE_level2_plot_df_ordered$level2_total_abundance_plot,
          probs = 0.95,
          na.rm = TRUE
        )
      ),
      oob = scales::squish,
      labels = scales::label_percent(
        accuracy = 0.01
      ),
      name = "Summed mean\nrelative abundance"
    ) +
    
    # Overlay bubbles.
    # Bubble size shows the mean abundance of the BRITE level 2 category.
    geom_point(
      aes(size = level2_mean_abundance_plot),
      shape = 21,
      color = "gray15",
      fill = "gray15",
      alpha = 0.65,
      stroke = 0.25
    ) +
    
    # Scale bubble area by mean relative abundance.
    # Values are shown as percentages to make small relative-abundance
    # values easier to interpret in the legend.
    scale_size_area(
      max_size = 4,
      limits = c(
        0,
        quantile(
          BRITE_level2_plot_df_ordered$level2_mean_abundance_plot,
          probs = 0.95,
          na.rm = TRUE
        )
      ),
      oob = scales::squish,
      breaks = c(
        0.0000,
        0.0001,
        0.00025,
        0.0005,
        0.0010,
        0.0020
      ),
      labels = scales::label_percent(
        accuracy = 0.01
      ),
      name = "Mean relative\nabundance"
    ) +
  
    # Split sample types by respiratory tract region
    facet_grid(
      KEGG_category_name ~ RT_category,
      scales = "free",
      space = "free",
      drop = TRUE,
      labeller = labeller(
        KEGG_category_name = label_wrap_gen(width = 14),
        RT_category = c(
          "Upper RT" = "Upper\nRT",
          "Intermediate RT" = "Intermediate\nRT",
          "Lower RT" = "Lower\nRT"
        )
      )
    ) +
    
    # Replace the internal unique y-axis IDs with clean KEGG subcategory labels.
    scale_y_discrete(
      labels = function(x) {
        sub(".*___", "", x)
      }
    ) +
    
    # Add publication-style labels
    labs(
      title = "Functional landscape across respiratory sample types",
      x = NULL,
      y = "KEGG BRITE Level 2 category"
    ) +
    
    # Use a clean base theme
    theme_bw() +
    
    # Refine axis text, title, grid, and legend appearance
    theme(
      
      # Rotate sample type labels for readability
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10,
        colour = "black"
      ),
      
      # Keep BRITE category labels compact
      axis.text.y = element_text(
        size = 10,
        colour = "black"
      ),
      
      # increase axis title size
      axis.title = element_text(size = 14, face = "bold"),
      
      # Make the title prominent but not oversized
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      
      # Remove major and minor grid lines behind tiles
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Keep legends on the right by default
      legend.position = "right",
      
      # Improve legend title readability
      legend.title = element_text(
        size = 14,
        face = "bold"
      ),
      # facet strip options
      strip.text.x = element_text(
        size = 8,
        face = "bold"
      ),
      
      strip.text.y.right = element_text(
        angle = 0,
        size = 11,
        face = "bold",
        margin = margin(
          t = 2,
          r = 4,
          b = 2,
          l = 4
        )
      ),
      strip.background = element_blank(),
      
      # Improve legend text readability
      legend.text = element_text(
        size = 12
      )
    )

# ---- 13.2 Save plot ----

  # save png
  ggsave("./Plots/Functional_landscape_BRITE_category_bubble_heatmap.png", 
         plot = BRITE_bubble_plot, 
         height = 24, 
         width = 24, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 600)
  
  # save svg
  ggsave("./Plots/Functional_landscape_BRITE_category_bubble_heatmap.svg", 
         plot = BRITE_bubble_plot, 
         height = 24, 
         width = 24, 
         units = "cm", 
         limitsize = FALSE)
  
  # save pdf
  ggsave(
  filename = "./Plots/Functional_landscape_BRITE_category_bubble_heatmap.pdf",
  plot = BRITE_bubble_plot,
  width = 9,
  height = 11
  )
  
  
# ---------------------------------------------------------------
# 12. Plot BRITE heatmap-bubble plot
# ---------------------------------------------------------------
  # Create subset of the above plot for main fig
  
  # ---- 12.2 Subset data ----
  
  KEGG_category_plot <- c("Metabolism",
                          "Genetic Information Processing",
                          "Environmental Information Processing"
                          )
  
  # subset plot df
  BRITE_level2_plot_df_sub <-
    BRITE_level2_plot_df_ordered %>%
    filter(KEGG_category_name %in% KEGG_category_plot)
  
### ---- 12.2 Generate ggplot ----
  
  # Initialize the plot with sample type on x-axis
  # and KEGG BRITE level 2 category on y-axis.
  
  BRITE_bubble_plot <- ggplot(
    BRITE_level2_plot_df_sub,
    aes(
      x = SampleTypev2,
      y = KEGG_subcategory_plot_id
    )
  )
  
  # build aesthetic details
  BRITE_bubble_plot <- 
    BRITE_bubble_plot +
    
    # Draw heatmap-style tiles.
    # Tile fill shows the total abundance of pathways/KOs
    # summarized into each BRITE level 2 category.
    geom_tile(
      aes(fill = level2_total_abundance_plot),
      color = "white",
      linewidth = 0.4
    ) +
    
    # Use a white-to-magenta gradient for Total mean_relative abundance.
    # The values are displayed as percentages because the original data
    # came from relative abundance profiles.
    scale_fill_gradient(
      low = "#FFF5FB",
      high = "#8B008B",
      trans = pseudo_log_trans(base = 10),
      limits = c(
        0,
        quantile(
          BRITE_level2_plot_df_sub$level2_total_abundance_plot,
          probs = 0.95,
          na.rm = TRUE
        )
      ),
      oob = scales::squish,
      labels = scales::label_percent(
        accuracy = 0.01
      ),
      name = "Cumulative mean\nrelative abundance"
    ) +
    
    # Overlay bubbles.
    # Bubble size shows the mean abundance of the BRITE level 2 category.
    geom_point(
      aes(size = level2_mean_abundance_plot),
      shape = 21,
      color = "gray15",
      fill = "gray15",
      alpha = 0.65,
      stroke = 0.25
    ) +
    
    # Scale bubble area by mean relative abundance.
    # Values are shown as percentages to make small relative-abundance
    # values easier to interpret in the legend.
    scale_size_area(
      max_size = 4,
      limits = c(
        0,
        quantile(
          BRITE_level2_plot_df_sub$level2_mean_abundance_plot,
          probs = 0.95,
          na.rm = TRUE
        )
      ),
      oob = scales::squish,
      breaks = c(
        0.0000,
        0.0001,
        0.00025,
        0.0005,
        0.0010,
        0.0020
      ),
      labels = scales::label_percent(
        accuracy = 0.01
      ),
      name = "Average mean\nrelative abundance"
    ) +
    
    # Split sample types by respiratory tract region
    facet_grid(
      KEGG_category_name ~ RT_category,
      scales = "free",
      space = "free",
      drop = TRUE,
      labeller = labeller(
        KEGG_category_name = label_wrap_gen(width = 14),
        RT_category = c(
          "Upper RT" = "Upper\nRT",
          "Intermediate RT" = "Intermediate\nRT",
          "Lower RT" = "Lower\nRT"
        )
      )
    ) +
    
    # Replace the internal unique y-axis IDs with clean KEGG subcategory labels.
    scale_y_discrete(
      labels = function(x) {
        sub(".*___", "", x)
      }
    ) +
    
    # Add publication-style labels
    labs(
      #title = "Functional landscape across respiratory sample types",
      x = NULL,
      y = "KEGG BRITE Level 2 category"
    ) +
    
    # ---- Legend controls ---- #
    guides(
      fill = guide_colourbar(
        order = 1,
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(6, "cm"),
        barheight = unit(0.4, "cm")
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        title.hjust = 0.5,
        nrow = 1,
        override.aes = list(
          shape = 21,
          fill = "gray15",
          colour = "gray15",
          alpha = 0.65,
          stroke = 0.25
        )
      )
    ) +
    
    # Use a clean base theme
    theme_bw() +
    
    # Refine axis text, title, grid, and legend appearance
    theme(
      
      # Rotate sample type labels for readability
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10,
        colour = "black"
      ),
      
      # Keep BRITE category labels compact
      axis.text.y = element_text(
        size = 10,
        colour = "black"
      ),
      
      # increase axis title size
      axis.title = element_text(size = 14, face = "bold"),
      
      # Make the title prominent but not oversized
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      
      # Remove major and minor grid lines behind tiles
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Keep legends on the bottom
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      
      # Improve spacing and alignment
      legend.box.just = "center",
      legend.spacing.x = unit(0.8, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
      
      # Improve legend title readability
      legend.title = element_text(
        size = 10,
        face = "bold"
      ),
      
      # Improve legend text readability
      legend.text = element_text(
        size = 8,
        face = "bold"
      ),
      
      # facet strip options
      strip.text.x = element_text(
        size = 8,
        face = "bold"
      ),
      
      strip.text.y.right = element_text(
        angle = 0,
        size = 11,
        face = "bold",
        margin = margin(
          t = 2,
          r = 4,
          b = 2,
          l = 4
        )
      ),
      strip.background = element_blank(),

      plot.margin = margin(t = 10,
                           l = 10,
                           r = 20
      )
      

    )
  
  # ---- 13.2 Save plot ----
  
  # save png
  ggsave("./Plots/Fig_3-D_Functional_landscape_BRITE_category_bubble_heatmap.png", 
         plot = BRITE_bubble_plot, 
         height = 16, 
         width = 20, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 600)
  
  # save svg
  ggsave("./Plots/Fig_3-D_Functional_landscape_BRITE_category_bubble_heatmap.svg", 
         plot = BRITE_bubble_plot, 
         height = 16, 
         width = 20, 
         units = "cm", 
         limitsize = FALSE)
  
  
  
  
  
  
  
  
  
  
  