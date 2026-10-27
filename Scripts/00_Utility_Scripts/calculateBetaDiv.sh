	#!/bin/bash 
printf "Program to calculate beta diversity from Metaphlan profiles.\n"
date

# load R module on FASRC
module load R

# init path variables
ProjectID='Merged_projects'
Date_Time=$(date +%d-%b-%H_%M)
WorkDir='/n/holyscratch01/huttenhower_lab/tshinde'
LogDir=${WorkDir}/Logs

# Give Metaphlan merged profile name her
Merged_profiles="${WorkDir}/Biobakery_run/Merged_profiles/Metaphlan4_taxonomic_profiles_merged.tsv"

# define which metrics you want to calculate
betaDiv_options=("bray-curtis" "jaccard" "clr" "aitchison")

# calculate them one by one
for betaDiv in "${betaDiv_options[@]}"; do
		printf '\nCalculating beta diversity metric "%s" for file %s\n' ${betaDiv} ${ProjectID}
		# call R script with metric param
		Rscript ${WorkDir}/Scripts/calculate_diversity.R -f ${Merged_profiles} -d beta -m ${betaDiv}

		printf "Results are stored in: diversity_analysis/*${betaDiv}.tsv \n"

		# sanity check
		ls -hs ${WorkDir}/diversity_analysis/*${betaDiv}.tsv
	done

printf "End program\n"
