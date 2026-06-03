# Code to merge Maaslin3 results form different Models into 1 file
# For LRT 
#    take: 'significant_results.tsv' from the "LRT_KO_Model-1" 
#    read and join (rowbind): all_results.tsv from each RT_KO_CF-Dataset-{N}_Model-2
#    output: "LRT_Maaslin3_merged_Coeff.txt"

# Load env
library(tidyverse)

# set wd
setwd("C:/Users/o_shinde/Desktop/MiPORT_Functions_DA/03_Humann/03_G_Maaslin3/LRT_KO_Models/")

# List of LRT datasets-Models 
Model_1_Name <- "LRT_KO_Model-1"
Model_2_Names <- list.files()[endsWith(list.files(), "Model-2")]

# sanity check
print("Joining results from the following models:")
print(Model_2_Names)

# -------------------------------
# 1. Read Model-1 significant results
# -------------------------------
M1_significant_Results <- read.table(
    paste0(Model_1_Name, "/significant_results.tsv"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE
)

# Assign metadata columns
M1_significant_Results$DatasetID <- "All"
M1_significant_Results$ModelID   <- "Model-1"
M1_significant_Results$SourceModel <- Model_1_Name

# -------------------------------
# 2. Read all Model-2 all_results.tsv files
# -------------------------------
Model2_res_List <- list()

for (eachM2 in Model_2_Names) {
    
    print(paste0("Reading file from model: ", eachM2))
    
    # get Dataset-ID value from folder name
    M2_Dataset_ID <- str_remove(eachM2,  "LRT_KO_") |> str_remove("_Model-2")
    
    # read file
    M2_Dataset_res <- read.table(
        paste0(eachM2, "/all_results.tsv"),
        sep = "\t",
        header = TRUE,
        check.names = FALSE
    )
    
    # Add metadata columns
    M2_Dataset_res$DatasetID   <- M2_Dataset_ID
    M2_Dataset_res$ModelID     <- "Model-2"
    M2_Dataset_res$SourceModel <- eachM2
    
    # Save in list
    Model2_res_List[[eachM2]] <- M2_Dataset_res
    
    print(dim(M2_Dataset_res))
}

# -------------------------------
# 3. Collapse all Model-2 results into one table
# -------------------------------
M2_all_results <- bind_rows(Model2_res_List)

print(dim(M2_all_results))

# sanity check
table(M2_all_results$DatasetID)

# -------------------------------
# 4. Standardize columns between Model-1 and Model-2
# -------------------------------
all_cols <- union(colnames(M1_significant_Results), colnames(M2_all_results))

M1_significant_Results <- M1_significant_Results %>%
    add_column(!!!setNames(
        rep(list(NA), length(setdiff(all_cols, colnames(.)))),
        setdiff(all_cols, colnames(.))
    ))

M2_all_results <- M2_all_results %>%
    add_column(!!!setNames(
        rep(list(NA), length(setdiff(all_cols, colnames(.)))),
        setdiff(all_cols, colnames(.))
    ))

# reorder columns identically
M1_significant_Results <- M1_significant_Results %>%
    select(all_of(all_cols))

M2_all_results <- M2_all_results %>%
    select(all_of(all_cols))

# -------------------------------
# 5. Merge Model-1 + Model-2 results
# -------------------------------
Merged_df <- bind_rows(
    M1_significant_Results,
    M2_all_results
)

# Move Model and Dataset column to the front
Merged_df <- Merged_df %>%
    relocate(all_of(c("DatasetID", "ModelID")))

print(dim(Merged_df))
print(table(Merged_df$ModelID, useNA = "ifany"))
print(table(Merged_df$DatasetID, useNA = "ifany"))

# -------------------------------
# 6. Optional: sort rows for readability
# -------------------------------
sort_cols <- intersect(c("feature", "metadata", "value", "DatasetID", "ModelID"), colnames(Merged_df))

if (length(sort_cols) > 0) {
    Merged_df <- Merged_df %>%
        arrange(across(all_of(sort_cols)))
}

# -------------------------------
# 7. Write merged output
# -------------------------------
write.table(
    Merged_df,
    file = "LRT_Maaslin3_merged_Coeff.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

print("Merged table written successfully:")
print("LRT_Maaslin3_merged_Coeff.txt")
