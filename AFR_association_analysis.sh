#!/bin/bash

# Navigate to the directory containing the filtered AFR_VCF file
cd ~/private/Project2/statSTR/AFR_calls_filtered

# Check if bcftools queries ran successfully
if [ $? -eq 0 ]; then

    # Summarize the genotype information in the filtered AFR_VCF file
    bcftools query -f "%CHROM\t%POS\t%PERIOD\t%END\t[%GB\t]\n" AFR_calls_filtered.vcf.gz > ~/private/Project2/Association_analysis/AFR_GB.txt

    # Summarize the sample names in the filtered AFR_VCF file
    bcftools query -l AFR_calls_filtered.vcf.gz > ~/private/Project2/Association_analysis/AFR_names.txt

    echo "bcftools queries completed successfully. Now running the association analysis script."

    # Navigate to the directory where the association analysis will be performed
    cd ~/private/Project2/Association_analysis/

    # Check if the python script ran successfully
    if [ $? -eq 0 ]; then
    	# Run the python script for association analysis
	python3 ~/public/project2-association/project2_association_script.py AFR 11 \
	~/private/Project2/Association_analysis/AFR_names.txt \
	~/private/Project2/Association_analysis/AFR_GB.txt > AFR_expr_results.tab

        echo "Association analysis complete. The results are saved in AFR_expr_results.tab."
    else
        echo "Association analysis script failed. Please check the input files and try again."
    fi
else
    echo "bcftools queries did not run successfully. Please check the AFR_VCF file and try again."
fi
