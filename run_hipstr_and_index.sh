#!/bin/bash
       
# Check if HipSTR ran successfully
if [ $? -eq 0 ]; then 
        # Run HipSTR for AFR_bam
    HipSTR --bam-files ~/private/Project2/Bam_file_list/AFR_bam_list/AFR_bam_files.txt \
           --fasta ~/public/genomes/Homo_sapiens_assembly38.fasta \
           --regions ~/private/Project2/hipstr_ref_small_bed/hipstr_ref_small.bed \
           --str-vcf ~/private/Project2/vcf/AFR_calls/AFR_calls.vcf.gz
           
         # Run HipSTR for EUR_bam
    HipSTR --bam-files ~/private/Project2/Bam_file_list/EUR_bam_list/EUR_bam_files.txt \
           --fasta ~/public/genomes/Homo_sapiens_assembly38.fasta \
           --regions ~/private/Project2/hipstr_ref_small_bed/hipstr_ref_small.bed \
           --str-vcf ~/private/Project2/vcf/EUR_calls/EUR_calls.vcf.gz
        
    echo "HipSTR ran successfully. Now indexing the VCF file."
   
    # Check if tabix ran successfully
    if [ $? -eq 0 ]; then
        # Index the output AFR VCF_call with tabix
        tabix -p vcf ~/private/Project2/vcf/AFR_calls/AFR_calls.vcf.gz

        # Index the output EUR VCF_call with tabix
        tabix -p vcf ~/private/Project2/vcf/EUR_calls/EUR_calls.vcf.gz
        
        echo "Indexing complete. The VCF file is now indexed."
    else
        echo "Indexing failed. Please check the VCF file and try again."
    fi
else
    echo "HipSTR did not run successfully. Please check the input files and parameters."
fi
