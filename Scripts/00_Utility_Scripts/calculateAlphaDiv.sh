#!/bin/bash 
printf "Program to calculate alpha diversity from Metaphlan profiles.\n"
date

# load R module on FASRC
module load R

# init path variables
ProjectID='Merged_projects'
Date_Time=$(date +%d-%b-%H_%M)
WorkDir='/n/holyscratch01/huttenhower_lab/tshinde'
LogDir=${WorkDir}/Logs

# Give Metaphlan merged profile name here
Merged_profiles="${WorkDir}/Metaphlan4_profiles/Merged_Profiles/Metaphlan4_taxonomic_profiles_merged.tsv"

# sanity check
ls -hs ${Merged_profiles}

# define which metrics you want to calculate
alphaDiv_options=("richness" "shannon" "simpson" "gini")

# calculate them one by one
for alphaDiv in "${alphaDiv_options[@]}"; do
		printf '\nCalculating alpha diversity metric "%s" for file %s\n' ${alphaDiv} ${ProjectID}
		
		# call R script with metric param
		Rscript ${WorkDir}/Scripts/calculate_diversity.R -f ${Merged_profiles} -d alpha -m ${alphaDiv}

		printf "Results are stored in: diversity_analysis/*${alphaDiv}.tsv \n"
		
		# sanity check
		ls -hs ${WorkDir}/diversity_analysis/*${alphaDiv}.tsv

	done

printf "End program\n"

