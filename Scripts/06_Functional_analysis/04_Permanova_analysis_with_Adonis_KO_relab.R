library(tidyverse)
library(vegan)

set.seed(8284)

# =========================
# 0. Setup
# =========================
OutDir <- "./03_E_Adonis"
dir.create(OutDir, showWarnings = FALSE)

nPerm <- 10000
nCores <- 10

# =========================
# 1. Read input
# =========================
KO_ext_relab <- read.table(
  "../03_BC_Humann4_merged_KO_extended_relab.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

Metadata_df <- read.table(
  "../../S3_MiPoRT_Metadata_global_QC_M4_RealAb_passed_samples_with_AN.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# =========================
# 2. Restrict to common samples
# =========================
Common_samples <- intersect(
  Metadata_df$SampleID,
  colnames(KO_ext_relab)[-1]
)

Metadata_Cols_of_Interest <- c(
  head(colnames(Metadata_df), 5),
  "SampleTypev2",
  "AgeGroup",
  "RT_category",
  "Sex"
)

Metadata_Cols_for_Adonis2 <- c(
  "BioProject",
  "SampleType",
  "Disease",
  "Healthy",
  "AgeGroup",
  "RT_category",
  "Sex"
)

Metadata_df <- Metadata_df %>%
  filter(SampleID %in% Common_samples) %>%
  select(all_of(Metadata_Cols_of_Interest)) %>%
  distinct(SampleID, .keep_all = TRUE)

KO_ext_relab <- KO_ext_relab %>%
  select(all_of(c("KO_ext_relab", Common_samples)))

# =========================
# 3. Recode metadata factors
# =========================
Metadata_df <- Metadata_df %>%
  mutate(
    RT_category = factor(
      RT_category,
      levels = c("URT", "IRT", "LRT"),
      labels = c(
        "Upper RT (N=1630)",
        "Intermediate RT (N=719)",
        "Lower RT (N=578)"
      )
    ),
    SampleTypev2 = factor(
      SampleTypev2,
      levels = c(
        "Nasal_Swab", "Buccal_mucosa", "Oral_swab", "Saliva",
        "Tongue_dorsum", "Supraglottal", "Sputum", "BAL", "Other"
      ),
      labels = c(
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
    )
  )

# Optional NA handling for categorical vars
Metadata_df <- Metadata_df %>%
  mutate(
    AgeGroup = replace_na(as.character(AgeGroup), "Unknown"),
    Sex      = replace_na(as.character(Sex), "Unknown"),
    AgeGroup = factor(AgeGroup),
    Sex      = factor(Sex, levels = c("Female", "Male", "Unknown"))
  )

# sanity check
table(Metadata_df$AgeGroup, useNA = 'ifany')
"
      Adult       Child Older_adult     Unknown Young_adult 
        269        1074         397         387         800
"
table(Metadata_df$Sex, useNA = 'ifany')
"
 Female    Male Unknown 
   1054    1367     506
"

# =========================
# 4. Build Bray-Curtis distance
# =========================
KO_ext_relab_t <- KO_ext_relab %>%
  column_to_rownames("KO_ext_relab") %>%
  t()

BC_dist <- vegdist(KO_ext_relab_t, method = "bray")

# =========================
# 5. Align metadata to distance matrix
# =========================
sample_ids <- attr(BC_dist, "Labels")

Metadata_df <- Metadata_df %>%
  filter(SampleID %in% sample_ids) %>%
  slice(match(sample_ids, SampleID))

rownames(Metadata_df) <- Metadata_df$SampleID

stopifnot(nrow(Metadata_df) == length(sample_ids))
stopifnot(all(Metadata_df$SampleID == sample_ids))

# =========================
# 6. Quick NA summary
# =========================
NA_summary <- tibble(
  Metadata = Metadata_Cols_for_Adonis2,
  N_available = sapply(
    Metadata_Cols_for_Adonis2,
    function(x) sum(!is.na(Metadata_df[[x]]))
  ),
  N_missing = sapply(
    Metadata_Cols_for_Adonis2,
    function(x) sum(is.na(Metadata_df[[x]]))
  )
)

print(NA_summary)

# =========================
# 7. Helper function for univariate PERMANOVA
# =========================
run_univariate_adonis <- function(var_name, dist_obj, meta_df, nPerm = 999, nCores = 8) {
  message("Running model: BC_dist ~ ", var_name)
  
  # use all rows currently available in metadata
  keep_idx <- !is.na(meta_df[[var_name]])
  
  meta_sub <- meta_df[keep_idx, , drop = FALSE]
  dist_sub <- as.dist(as.matrix(dist_obj)[keep_idx, keep_idx])
  
  # skip impossible models
  if (length(unique(meta_sub[[var_name]])) < 2) {
    return(tibble(
      Metadata = var_name,
      N = nrow(meta_sub),
      R2 = NA_real_,
      F = NA_real_,
      p_value = NA_real_,
      Significance = NA_character_
    ))
  }
  
  model_formula <- as.formula(paste("dist_sub ~", var_name))
  
  fit <- adonis2(
    formula = model_formula,
    data = meta_sub,
    permutations = nPerm,
    parallel = nCores
  )
  
  tibble(
    Metadata = var_name,
    N = nrow(meta_sub),
    R2 = fit$R2[1],
    F = fit$F[1],
    p_value = fit$`Pr(>F)`[1],
    Significance = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE             ~ "NS"
    )
  )
}

# =========================
# 8. Run all univariate models
# =========================
Adonis_results <- map_dfr(
  Metadata_Cols_for_Adonis2,
  ~ run_univariate_adonis(
    var_name = .x,
    dist_obj = BC_dist,
    meta_df = Metadata_df,
    nPerm = nPerm,
    nCores = nCores
  )
) %>%
  mutate(Statistic_type = "Univariate") %>%
  arrange(desc(R2))

print(Adonis_results)
Adonis_results <- Adonis_results %>% filter(Metadata != 'PatientID')
 
# =========================
# 9. Multi-variate PERMANOVA Tests
# =========================

# Permanova with all variables in Metadata_Cols_for_Adonis2 
adonis2_multivar_all <- 
  adonis2(BC_dist ~ 
            SampleType + BioProject + Disease + 
            RT_category + AgeGroup + Healthy + Sex, 
          data = Metadata_df[,Metadata_Cols_of_Interest],
          permutations = nPerm, 
          parallel = 8)

adonis2_multivar_all

# save this multi variable result plot
Multivar_model <- tibble(
  Metadata = 'All',
  N = 2927,
  R2 = adonis2_multivar_all$R2[1],
  F = adonis2_multivar_all$F[1],
  p_value = adonis2_multivar_all$`Pr(>F)`[1],
  Significance = case_when(
    p_value <= 0.001 ~ "***",
    p_value <= 0.01  ~ "**",
    p_value <= 0.05  ~ "*",
    TRUE             ~ "NS"
  ),
  Statistic_type = "Multivariate"
  )

# Add this model result to Univar result df
Adonis_results_merged <- rbind(Multivar_model, Adonis_results)

# sanity check
str(Adonis_results_merged)

# =========================
# 10. Save output
# =========================
write.table(
  Adonis_results_merged,
  file = file.path(OutDir, "KO_univariate_adonis_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# Plot: Barplot Rsq with p-values as text annotations
  ggPlotObj <- ggplot(Adonis_results_merged, 
                      aes(x = reorder(Metadata, -R2), 
                          y = R2*100)) 
  
  # save plot
  savePlot <- ggPlotObj +
    geom_col(fill = "gray20", 
             alpha = 0.8,
             width = 0.6) +  # Barplot for Rsq
    geom_text(aes(label = Significance), # Add significance stars
              hjust = 0.5, vjust = 0.2,
              size = 6, 
              color = "gray20") + 
    geom_text(aes(label = round(R2*100, 1)), # Add significance stars
              hjust = 0.5, vjust = 1.5,
              size = 9, 
              color = "white") +
    labs(
      title = "Functional variance explained by metadata",
      subtitle = "KO-level",
      x = "Metadata",
      y = "R-squared (%)",
      caption = paste0(paste0("n permutations = ", nPerm)," || Significance stars: ***p ≤ 0.001, **p ≤ 0.01, *p ≤ 0.05")
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 18, face = "bold", colour = "Black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 20, face = "bold"),
      title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 18),
      plot.caption = element_text(size = 10)
    )
  
  savePlot
  
  # save png plot
  ggsave("./03_E_Adonis/All_model_PERMANOVA_KO_Batch_corrected.png",
         plot= savePlot, 
         height = 26, 
         width = 30, 
         units = "cm", 
         limitsize = FALSE, 
         dpi = 300)
  
  # save svg plot
  ggsave("./03_E_Adonis/All_model_PERMANOVA_KO_Batch_corrected.svg",
         plot=savePlot, 
         height = 18, 
         width = 28, 
         units = "cm", 
         limitsize = FALSE)
  
  
  