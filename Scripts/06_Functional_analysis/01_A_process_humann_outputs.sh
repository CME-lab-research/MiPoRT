#!/usr/bin/env bash

# ============================================================
# 1. Script overview
# ============================================================
# This script merges and pre-processes HUMAnN-derived functional profiles
# generated for the MiPoRT respiratory metagenome atlas.
#
# Main outputs:
#   1. Merged HUMAnN gene-family table
#   2. Gene-family table without HUMAnN bucket features
#   3. Unstratified gene-family table
#   4. KEGG Ortholog (KO) tables using standard and expanded UniRef90-to-KO mappings
#   5. eggNOG-regrouped table
#   6. Relative-abundance-normalized KO tables
#   7. Merged and normalized MetaCyc pathway-abundance table
#   8. Per-sample KO mapping coverage summaries
#
# Notes for reuse:
#   - Replace the path variables below before running this script.
#   - This script assumes one HUMAnN output directory per BioProject/cohort.
#   - Large memory is required for HUMAnN regrouping at MiPoRT scale.
#   - Some zero-signal samples were removed using separate R scripts before renormalization.

set -euo pipefail

printf '\n---------------------- Start program ----------------------\n'
echo "$(date +%d-%b\|%H:%M)"

# ============================================================
# 2. Configure paths, environment, and cohort list
# ============================================================

# ---- 2.1 Define conda/mamba environment variables ----
# Define the directory containing the mamba/conda activation script.
# This is kept as a user-editable path for portability across HPC systems.
MambaDir="/path/to/conda_mamba/bin"

# Define the conda environment that contains HUMAnN and MetaPhlAn.
CondaEnvName="BBHm_v4"

# ---- 2.2 Define project working directories ----
# Define the main working directory containing all HUMAnN output folders.
WorkDir="/path/to/MiPoRT_functional_analysis"

# Define the directory where copied gene-family profiles will be pooled.
gFam_dataDir="${WorkDir}/Pooled_geneFamily_profiles"

# Define the directory where copied pathway profiles will be pooled.
Pathway_dataDir="${WorkDir}/Pooled_pathway_profiles"

# Define the HUMAnN utility-mapping directory.
# This is needed only for the expanded UniRef90-to-KO mapping file.
UtilDB="/path/to/humann/utility_mapping"

# ---- 2.3 Define cohort-specific HUMAnN output folders ----
# These folders contain per-sample HUMAnN output files.
# Expected file patterns:
#   - *_2_genefamilies.tsv
#   - *_4_pathabundance.tsv
HumannProfileDirs=(
  "LRT_Mix_B7_Humann_profiles"
  "PRJEB27079_Humann_profiles"
  "PRJEB29011_Humann_profiles"
  "PRJNA413615_Humann_profiles"
  "PRJNA470402_Humann_profiles"
  "PRJNA48479_Humann_profiles"
  "PRJNA494034_Humann_profiles"
  "PRJNA687506_Humann_profiles"
  "PRJNA756530_Humann_profiles"
  "PRJNA757846_Humann_profiles"
  "PRJNA762218_Humann_profiles"
  "PRJNA917836_Humann_profiles"
  "UniBergen_Humann_profiles"
)

# ---- 2.4 Move to the working directory and create output directories ----
# Move into the working directory so all relative paths resolve correctly.
cd "${WorkDir}"

# Create pooled input directories if they do not already exist.
mkdir -vp "${gFam_dataDir}" "${Pathway_dataDir}"

# ============================================================
# 3. Activate software environment and print tool versions
# ============================================================

# ---- 3.1 Activate the HUMAnN analysis environment ----
# Activate the environment used during the MiPoRT HUMAnN analysis.
source "${MambaDir}/activate" "${CondaEnvName}"

# ---- 3.2 Print software versions for reproducibility ----
# Print versions without stopping the script if version commands are unavailable.
printf "\nEnvironment details:\n"
metaphlan --version || true
humann --version || true

# ============================================================
# 4. Merge HUMAnN gene-family profiles
# ============================================================

# ---- 4.1 Count input gene-family profiles ----
# Count all gene-family profiles detected across cohort-specific HUMAnN folders.
printf "\nCounting HUMAnN gene-family profiles...\n"
ls ./*_Humann_profiles/*2_genefamilies.tsv | wc -l

# ---- 4.2 Record cohort-level input counts ----
# These checks document expected file counts used in the MiPoRT revision analysis.
# Expected counts from the original run:
#   LRT_Mix_B7: 316
#   PRJEB27079: 4
#   PRJEB29011: 204
#   PRJNA413615: 271
#   PRJNA470402: 58
#   PRJNA48479: 1116
#   PRJNA494034: 82
#   PRJNA687506: 270
#   PRJNA756530: 33
#   PRJNA757846: 113
#   PRJNA762218: 23
#   PRJNA917836: 760
#   UniBergen: 253
for profile_dir in "${HumannProfileDirs[@]}"; do
  printf "%s\t" "${profile_dir}"
  ls "${profile_dir}"/*2_genefamilies.tsv | wc -l
done

# ---- 4.3 Copy gene-family profiles into one pooled directory ----
# Copy all gene-family profiles into one directory before using humann_join_tables.
printf "\nCopying gene-family profiles to: %s\n" "${gFam_dataDir}"
for profile_dir in "${HumannProfileDirs[@]}"; do
  rsync -a -P "${profile_dir}"/*2_genefamilies.tsv "${gFam_dataDir}/"
done

# ---- 4.4 Check pooled gene-family profile count ----
# Expected original MiPoRT count: 3503 samples.
printf "\nPooled gene-family profile count:\n"
ls "${gFam_dataDir}"/*2_genefamilies.tsv | wc -l

# ---- 4.5 Join gene-family profiles into one feature table ----
# Merge all per-sample HUMAnN gene-family profiles into a single wide table.
humann_join_tables \
  --verbose \
  --input "${gFam_dataDir}/" \
  --output Humann4_merged_genefamilies.tsv \
  --file_name 2_genefamilies.tsv

# ---- 4.6 Check merged gene-family table size ----
# Expected original MiPoRT output:
#   File size: ~17 GB
#   Row count: 12,050,103 rows
ls -hs Humann4_merged_genefamilies.tsv
wc -l Humann4_merged_genefamilies.tsv

# ============================================================
# 5. Clean and unstratify gene-family profiles
# ============================================================

# ---- 5.1 Remove HUMAnN bucket features ----
# Remove non-biological HUMAnN summary rows before regrouping.
# Examples: READS_UNMAPPED and UniRef50_unknown.
grep -v -E '^READS_UNMAPPED\b|^UniRef50_unknown\b' \
  Humann4_merged_genefamilies.tsv \
  > Humann4_merged_genefamilies_no_buckets.tsv

# Expected original MiPoRT row count: 12,050,102 rows.
wc -l Humann4_merged_genefamilies_no_buckets.tsv

# ---- 5.2 Keep only unstratified gene-family rows ----
# HUMAnN stratified rows contain a pipe character followed by taxonomic labels.
# These rows are removed so each feature represents the total abundance across taxa.
awk -F'\t' 'NR==1 || $1 !~ /\|/' \
  Humann4_merged_genefamilies_no_buckets.tsv \
  > Humann4_merged_genefamilies_no_buckets_unstrat.tsv

# Expected original MiPoRT row count: 5,364,641 rows.
wc -l Humann4_merged_genefamilies_no_buckets_unstrat.tsv

# ============================================================
# 6. Regroup UniRef90 gene families to functional annotations
# ============================================================

# ---- 6.1 Regroup UniRef90 features to KEGG Orthologs using the standard mapping ----
# This creates the standard HUMAnN KO table.
# Original MiPoRT resource note: ~400 GB resident memory was required.
humann_regroup_table \
  --input Humann4_merged_genefamilies_no_buckets_unstrat.tsv \
  --groups uniref90_ko \
  --output Humann4_merged_KO.tsv

# Original MiPoRT regrouping summary:
#   Original feature count: 5,364,640
#   Grouped 1+ times: 231,954 (4.3%)
#   Grouped 2+ times: 1,191 (0.0%)

# ---- 6.2 Regroup UniRef90 features to KEGG Orthologs using the expanded mapping ----
# This creates the expanded KO table used to improve KO recovery.
# Original MiPoRT resource note: ~600 GB resident memory was required.
humann_regroup_table \
  --input Humann4_merged_genefamilies_no_buckets_unstrat.tsv \
  --custom "${UtilDB}/map_ko-expanded_uniref90.txt.gz" \
  --output Humann4_merged_KO_expanded.tsv

# Original MiPoRT regrouping summary:
#   Original feature count: 5,364,640
#   Grouped 1+ times: 951,797 (17.7%)
#   Grouped 2+ times: 34,829 (0.6%)

# ---- 6.3 Regroup UniRef90 features to eggNOG annotations ----
# This table was generated as an additional functional annotation layer.
# Original MiPoRT resource note: ~600 GB resident memory was required.
humann_regroup_table \
  --input Humann4_merged_genefamilies_no_buckets_unstrat.tsv \
  --groups uniref90_eggnog \
  --output Humann4_merged_eggNOG.tsv

# Original MiPoRT regrouping summary:
#   Original feature count: 5,364,640
#   Grouped 1+ times: 328,366 (6.1%)
#   Grouped 2+ times: 313,493 (5.8%)

# ---- 6.4 Check regrouped table dimensions ----
# Expected original MiPoRT row counts:
#   Humann4_merged_genefamilies_no_buckets_unstrat.tsv: 5,364,641
#   Humann4_merged_KO.tsv: 10,665
#   Humann4_merged_KO_expanded.tsv: 11,951
#   Humann4_merged_eggNOG.tsv: 74,729
wc -l Humann4_merged_genefamilies_no_buckets_unstrat.tsv
wc -l Humann4_merged_KO.tsv
wc -l Humann4_merged_KO_expanded.tsv
wc -l Humann4_merged_eggNOG.tsv

# ============================================================
# 7. Filter zero-signal KO samples and renormalize KO tables
# ============================================================

# ---- 7.1 Document manual/intermediate zero-sample filtering ----
# Some samples had zero KO signal after regrouping.
# These were removed before relative-abundance normalization.
# In the original analysis, this filtering was done using separate R scripts:
#   - 03_B_filter_zero_KO_samples.R
#
# Expected filtered input files:
#   - 01_Humann4_merged_KO_filtered.tsv
#   - 01_Humann4_merged_KO_extended_filtered.tsv

# ---- 7.2 Renormalize the standard KO table to relative abundance ----
# Renormalization prepares the table for downstream MMUPHin/MaAsLin3 analyses.
humann_renorm_table \
  --input 01_Humann4_merged_KO_filtered.tsv \
  --output 02_Humann4_merged_KO_relab.tsv \
  --units relab

# ---- 7.3 Check column consistency before and after standard KO renormalization ----
# These checks confirm that the number of sample columns is preserved.
head -1 01_Humann4_merged_KO_filtered.tsv | tr -cd '\t' | wc -c | awk '{print $1+1}'
head -2 01_Humann4_merged_KO_filtered.tsv | tail -n 1 | tr -cd '\t' | wc -c | awk '{print $1+1}'
head -1 02_Humann4_merged_KO_relab.tsv | tr -cd '\t' | wc -c | awk '{print $1+1}'
head -2 02_Humann4_merged_KO_relab.tsv | tail -n 1 | tr -cd '\t' | wc -c | awk '{print $1+1}'

# ---- 7.4 Renormalize the expanded KO table to relative abundance ----
# The expanded KO table uses the custom mapping and was the preferred KO table for revision analyses.
humann_renorm_table \
  --input 01_Humann4_merged_KO_extended_filtered.tsv \
  --output 02_Humann4_merged_KO_extended_relab.tsv \
  --units relab

# ---- 7.5 Check column consistency before and after expanded KO renormalization ----
# Expected original MiPoRT column count after zero-sample filtering: 3136 columns.
head -1 01_Humann4_merged_KO_extended_filtered.tsv | tr -cd '\t' | wc -c | awk '{print $1+1}'
head -2 01_Humann4_merged_KO_extended_filtered.tsv | tail -n 1 | tr -cd '\t' | wc -c | awk '{print $1+1}'
head -1 02_Humann4_merged_KO_extended_relab.tsv | tr -cd '\t' | wc -c | awk '{print $1+1}'
head -2 02_Humann4_merged_KO_extended_relab.tsv | tail -n 1 | tr -cd '\t' | wc -c | awk '{print $1+1}'

# ============================================================
# 8. Estimate per-sample KO mapping coverage
# ============================================================

# ---- 8.1 Calculate total gene-family abundance per sample ----
# This provides the denominator for KO mapping coverage.
awk -F'\t' '
NR==1 {
  for (i=2; i<=NF; i++) name[i]=$i
  next
}
{
  for (i=2; i<=NF; i++) s[i]+=$i
}
END {
  for (i=2; i<=NF; i++) print name[i], s[i]
}' Humann4_merged_genefamilies_no_buckets_unstrat.tsv \
  > Total_sample_gene_abundance.txt

# ---- 8.2 Calculate total standard KO abundance per sample ----
# This provides the standard KO numerator for each sample.
awk -F'\t' '
NR==1 {
  for (i=2; i<=NF; i++) name[i]=$i
  next
}
{
  for (i=2; i<=NF; i++) s[i]+=$i
}
END {
  for (i=2; i<=NF; i++) print name[i], s[i]
}' Humann4_merged_KO.tsv \
  > Total_sample_KO_abundance.txt

# ---- 8.3 Calculate total expanded KO abundance per sample ----
# This provides the expanded KO numerator for each sample.
awk -F'\t' '
NR==1 {
  for (i=2; i<=NF; i++) name[i]=$i
  next
}
{
  for (i=2; i<=NF; i++) s[i]+=$i
}
END {
  for (i=2; i<=NF; i++) print name[i], s[i]
}' Humann4_merged_KO_expanded.tsv \
  > Total_sample_KO_abundance_with_ExtendedMapping.txt

# ---- 8.4 Compute KO mapping coverage percentages ----
# This compares mapped KO abundance against total gene-family abundance.
# Typical KO coverage in this analysis was expected to be approximately 60-70%.
join <(sort Total_sample_gene_abundance.txt) \
  <(sort Total_sample_KO_abundance.txt) \
  | join - <(sort Total_sample_KO_abundance_with_ExtendedMapping.txt) \
  | awk 'BEGIN {OFS="\t"; print "SampleID", "Standard_KO_percent", "Expanded_KO_percent"}
{
  sample=$1
  gene=$2
  ko_std=$3
  ko_ext=$4

  pct_std=(ko_std/gene)*100
  pct_ext=(ko_ext/gene)*100

  printf "%s\t%.2f\t%.2f\n", sample, pct_std, pct_ext
}' \
  > KO_mapping_coverage_percent_per_sample.tsv

# ============================================================
# 9. Merge HUMAnN pathway-abundance profiles
# ============================================================

# ---- 9.1 Count input pathway-abundance profiles ----
# Count all pathway-abundance profiles detected across cohort-specific HUMAnN folders.
printf "\nCounting HUMAnN pathway-abundance profiles...\n"
ls ./*_Humann_profiles/*4_pathabundance.tsv | wc -l

# ---- 9.2 Record cohort-level pathway profile counts ----
# These checks document expected file counts used in the MiPoRT revision analysis.
# Expected counts match the gene-family profile counts above.
for profile_dir in "${HumannProfileDirs[@]}"; do
  printf "%s\t" "${profile_dir}"
  ls "${profile_dir}"/*4_pathabundance.tsv | wc -l
done

# ---- 9.3 Check pathway feature counts per sample ----
# Some samples can have very low or zero pathway signal.
# Original MiPoRT note: 183 samples had almost no pathway counts.
wc -l \
  LRT_Mix_B7_Humann_profiles/*4_pathabundance.tsv \
  PRJEB27079_Humann_profiles/*4_pathabundance.tsv \
  PRJEB29011_Humann_profiles/*4_pathabundance.tsv \
  PRJNA413615_Humann_profiles/*4_pathabundance.tsv \
  PRJNA470402_Humann_profiles/*4_pathabundance.tsv \
  PRJNA48479_Humann_profiles/*4_pathabundance.tsv \
  PRJNA494034_Humann_profiles/*4_pathabundance.tsv \
  PRJNA687506_Humann_profiles/*4_pathabundance.tsv \
  PRJNA756530_Humann_profiles/*4_pathabundance.tsv \
  PRJNA757846_Humann_profiles/*4_pathabundance.tsv \
  PRJNA762218_Humann_profiles/*4_pathabundance.tsv \
  PRJNA917836_Humann_profiles/*4_pathabundance.tsv \
  UniBergen_Humann_profiles/*4_pathabundance.tsv \
  > Pathway_feature_counts_SampleWise.txt

# ---- 9.4 Copy pathway-abundance profiles into one pooled directory ----
# Copy all pathway-abundance profiles into one directory before using humann_join_tables.
printf "\nCopying pathway-abundance profiles to: %s\n" "${Pathway_dataDir}"
for profile_dir in "${HumannProfileDirs[@]}"; do
  rsync -a -P "${profile_dir}"/*4_pathabundance.tsv "${Pathway_dataDir}/"
done

# ---- 9.5 Check pooled pathway profile count ----
# Expected original MiPoRT count: 3503 samples.
printf "\nPooled pathway-abundance profile count:\n"
ls "${Pathway_dataDir}"/*4_pathabundance.tsv | wc -l

# ---- 9.6 Join pathway-abundance profiles into one feature table ----
# Merge all per-sample HUMAnN pathway-abundance profiles into a single wide table.
humann_join_tables \
  --verbose \
  --input "${Pathway_dataDir}/" \
  --output Humann4_merged_pathways.tsv \
  --file_name 4_pathabundance.tsv

# Expected original MiPoRT row count: 79,044 rows.
wc -l Humann4_merged_pathways.tsv

# ============================================================
# 10. Clean, filter, and renormalize pathway-abundance profiles
# ============================================================

# ---- 10.1 Remove pathway bucket features and stratified rows ----
# Remove HUMAnN pathway summary rows and taxon-stratified rows.
awk -F'\t' '
NR==1 {print; next}
$1 !~ /^UNMAPPED/ &&
$1 !~ /^UNINTEGRATED/ &&
$1 !~ /\|/
' Humann4_merged_pathways.tsv \
  > Humann4_merged_pathways_no_buckets.tsv

# ---- 10.2 Confirm pathway bucket and stratified rows were removed ----
# These commands should return no matching rows after filtering.
grep -E '^UNMAPPED|^UNINTEGRATED' Humann4_merged_pathways_no_buckets.tsv || true
grep '\|' Humann4_merged_pathways_no_buckets.tsv || true

# ---- 10.3 Check pathway table dimensions and header ----
# Expected original MiPoRT row count after filtering: 806 rows, including header.
head -n 1 Humann4_merged_pathways_no_buckets.tsv
wc -l Humann4_merged_pathways_no_buckets.tsv

# ---- 10.4 Check for negative pathway abundance values ----
# Negative values are not expected in HUMAnN pathway-abundance outputs.
awk -F'\t' '
NR>1 {
  for (i=2; i<=NF; i++) if ($i<0) print "Negative:", $1
}' Humann4_merged_pathways_no_buckets.tsv

# ---- 10.5 Count non-zero pathway features per sample ----
# This identifies samples with little or no pathway signal.
awk -F'\t' '
NR>1 {
  for (i=2; i<=NF; i++) if ($i>0) count[i]++
}
END {
  for (i=2; i<=NF; i++) print "Column", i, "non-zero rows:", count[i]
}' Humann4_merged_pathways_no_buckets.tsv \
  > Path_empty_samples.txt

# ---- 10.6 Document manual/intermediate zero-sample filtering for pathways ----
# Some samples had all-zero pathway features.
# These samples were removed before relative-abundance normalization.
# In the original analysis, this filtering was done using a separate R script:
#   - 03_filter_zero_PathAb_samples.R
#
# Expected filtered input file:
#   - Humann4_merged_pathways_no_buckets_no_ZeroSamples.tsv

# ---- 10.7 Renormalize the filtered pathway table to relative abundance ----
# The --special n flag prevents HUMAnN special features from being retained.
# The --update-snames flag standardizes sample names in the output table.
humann_renorm_table \
  --input Humann4_merged_pathways_no_buckets_no_ZeroSamples.tsv \
  --output Humann4_merged_pathways_relab.tsv \
  --special n \
  --units relab \
  --update-snames

# ============================================================
# 11. Close software environment
# ============================================================

# ---- 11.1 Deactivate conda environment ----
# Deactivate the environment without failing if conda is unavailable.
conda deactivate || true

# ---- 11.2 Print completion timestamp ----
printf '\n----------------------  End program  ----------------------\n'
echo "$(date +%d-%b\|%H:%M)"
