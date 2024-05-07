# Read the CSV file which is the taxonomy file (genomes-all_metadata) downloaded from UHGG, assuming it is comma-separated
metadata_taxonomy <- read.csv("metadata_taxonomy_of_the_genome_from_UHGG_v.2.0.2.csv", header = TRUE)

# Check if Lineage ends with ";g__;s__" and has nothing after it
to_modify_g_s <- grepl(";g__;s__$", metadata_taxonomy$Lineage)

# Loop through the rows and add Species_rep after "g__" and ";s__" in the Lineage for rows that meet the condition
for (i in which(to_modify_g_s)) {
  metadata_taxonomy$Lineage[i] <- gsub(";g__", paste0(";g__", metadata_taxonomy$Species_rep[i]), metadata_taxonomy$Lineage[i])
  metadata_taxonomy$Lineage[i] <- gsub(";s__", paste0(";s__", metadata_taxonomy$Species_rep[i]), metadata_taxonomy$Lineage[i])
}


# Check if Lineage ends with "s__" and has nothing after it
to_modify <- grepl("^.*;s__$", metadata_taxonomy$Lineage)

# Add Species_rep to the end of Lineage for rows that meet the condition
metadata_taxonomy$Lineage[to_modify] <- paste0(metadata_taxonomy$Lineage[to_modify], metadata_taxonomy$Species_rep[to_modify])

# Print the updated data
print(metadata_taxonomy)
#write.csv(metadata_taxonomy, "metadata_taxonomy.csv", row.names = FALSE)

# Replace spaces with underscores in Lineage column
metadata_taxonomy$Lineage <- gsub(" ", "_", metadata_taxonomy$Lineage)


# Read the CSV file, assuming it is comma-separated
filtered_data <- read.csv("filtered_data_species_level.csv", header = TRUE)

# Replace "k__" with "d__" and "|" with ";" in the first column (OTU_ID)
filtered_data$OTU_ID <- gsub("k__", "d__", filtered_data$OTU_ID)
filtered_data$OTU_ID <- gsub("\\|", ";", filtered_data$OTU_ID)

# Print the updated data
print(filtered_data)


# Match OTU_ID values between the two data frames
match_indices <- match(filtered_data$OTU_ID, metadata_taxonomy$Lineage)

# Retrieve the Length values based on the matched indices
filtered_data$Length <- metadata_taxonomy$Length[match_indices]

# Add the Length column as the second column in filtered_data
filtered_data <- cbind(filtered_data[,1], Length = filtered_data$Length, filtered_data[,2:ncol(filtered_data)])

# Print the updated data
print(filtered_data)

#################Calculation of RPKM######################
# Calculate total reads per sample
total_reads <- colSums(filtered_data[, -c(1, 2)])

# Calculate per million scaling factor (sequencing depth)
scaling_factor <- total_reads / 1e6

# Normalize read counts by the scaling factor
normalized_data <- filtered_data[, -c(1, 2)] / scaling_factor

# Calculate RPKM values (based on Genome Length, kilobases)
rpkm_values <- normalized_data / (filtered_data$Length / 1000)

# Add the gene and length columns back
rpkm_data <- cbind(filtered_data[, c(1, 2)], rpkm_values)

# Write the result to a new CSV file
write.csv(rpkm_data, "rpkm_output.csv", row.names = FALSE)
