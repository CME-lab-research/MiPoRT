# Script for a basic workflow to do QC of paired-end (PE) sequencing reads for microbiome wgs/amplicon samples. 
# 22/03/2024
# Tejus Shinde (MEL)

# QC process pipeline
	# check quality before QC
	1. fastqc -t=2 ${R1_File} ${R2_File} -o ${OutDir}

	# remove host-DNA reads with BBTools
	2.  removehuman.sh in=${R1_File} outu=${R1_out} t=${CPUs_per_job}
	2.1 removehuman.sh in=${R2_File} outu=${R2_out} t=${CPUs_per_job}

	# Correct bases, remove and filter bad-quality reads
	3. fastp --in1 ${R1_in_fastqName} --in2 ${R2_in_fastqName} --out1 ${R1_out_fastqName} --out2 ${R2_out_fastqName} --average_qual 20 --detect_adapter_for_pe --correction --overrepresentation_analysis --low_complexity_filter 

	# Repair mate-pairs
	4. bash ${BBMAP_dir}/repair.sh in1=${R1_out_fastqName} in2=${R2_out_fastqName} out1=${R1_QC_File} out2=${R2_QC_File} repair ignorebadquality=t tossbrokenreads=t overwrite=t -Xmx1G

	# check quality before QC
	5. fastqc -t=2 ${R1_QC_File} ${R2_QC_File} -o ${OutDir}

	# pool fastqc reports with multiqc
	6. multiqc fastqc_PRJxx_QC/*_{1,2}_fastqc.zip --outdir="multiqc_${ProjectID}_QC"/

# For efficiency and multi sample analysis these cmds are added to cluster specific scripts and run.
