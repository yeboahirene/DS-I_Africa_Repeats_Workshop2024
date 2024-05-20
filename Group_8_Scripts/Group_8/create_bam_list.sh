#!/bin/bash

# Create a text file with a list of AFR_BAMs to genotype
ls ~/public/project2-association/dataset/AFR_bam/*bam > \
~/private/Project2/Bam_file_list/AFR_bam_list/AFR_bam_files.txt

# Create a text file with a list of EUR_BAMs to genotype
ls ~/public/project2-association/dataset/EUR_bam/*bam > \
~/private/Project2/Bam_file_list/EUR_bam_list/EUR_bam_files.txt

# Print a message to indicate completion
echo "List of BAM files created successfully."
