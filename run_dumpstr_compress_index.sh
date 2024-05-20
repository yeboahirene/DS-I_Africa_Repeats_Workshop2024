#!/bin/bash

# Check if dumpSTR ran successfully
if [ $? -eq 0 ]; then

    # Run dumpSTR with the specified options
    dumpSTR --vcf ~/private/Project2/vcf/AFR_calls/AFR_calls.vcf.gz \
            --out ~/private/Project2/statSTR/AFR_calls_filtered/AFR_calls_filtered \
            --vcftype hipstr \
            --hipstr-min-call-Q 0.9 \
            --hipstr-min-call-DP 10
            
    # Run dumpSTR with the specified options
    dumpSTR --vcf ~/private/Project2/vcf/EUR_calls/EUR_calls.vcf.gz \
            --out ~/private/Project2/statSTR/EUR_calls_filtered/EUR_calls_filtered \
            --vcftype hipstr \
            --hipstr-min-call-Q 0.9 \
            --hipstr-min-call-DP 10
            
   echo "dumpSTR ran successfully. Now compressing the VCF file."

    # Check if bgzip ran successfully
    if [ $? -eq 0 ]; then
    
        # Compress the filtered AFR_VCF file with bgzip
        bgzip ~/private/Project2/statSTR/AFR_calls_filtered/AFR_calls_filtered.vcf

        # Compress the filtered EUR_VCF file with bgzip
        bgzip ~/private/Project2/statSTR/EUR_calls_filtered/EUR_calls_filtered.vcf

        echo "Compression complete. Now indexing the VCF file."

        # Check if tabix ran successfully
        if [ $? -eq 0 ]; then
        
            # Index the compressed AFR_VC file with tabix
            tabix -p vcf ~/private/Project2/statSTR/AFR_calls_filtered/AFR_calls_filtered.vcf.gz

            # Index the compressed EUR_VCF file with tabix
            tabix -p vcf ~/private/Project2/statSTR/EUR_calls_filtered/EUR_calls_filtered.vcf.gz
            
            # Inspect repeat length distributions by population for AFR_filtered_vcf_calls  using statSTR
            statSTR --vcf ~/private/Project2/statSTR/AFR_calls_filtered/AFR_calls_filtered.vcf.gz \
            --plot-afreq \
            --out ~/private/Project2/statSTR/AFR_calls_filtered/AFR
            
            # Inspect repeat length distributions by population for EUR_filtered_vcf_calls  using statSTR
            statSTR --vcf ~/private/Project2/statSTR/EUR_calls_filtered/EUR_calls_filtered.vcf.gz \
            --plot-afreq \
            --out ~/private/Project2/statSTR/EUR_calls_filtered/EUR
        
            echo "Indexing complete. The VCF file is now indexed."
        else
            echo "Indexing failed. Please check the compressed VCF file and try again."
        fi
    else
        echo "Compression failed. Please check the VCF file and try again."
    fi
else
    echo "dumpSTR did not run successfully. Please check the input VCF file and parameters."
fi
