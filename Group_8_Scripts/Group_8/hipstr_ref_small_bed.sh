#!/bin/bash

# Filter STRs within chr11:57520464-57529754
awk '($1=="chr11" && $2 >=57520464 && $3 <= 57529754)' ~/public/project2-association/dataset/hg38.hipstr_reference.bed > ~/private/Project2/hipstr_ref_small_bed/hipstr_ref_small.bed

# Print a message to indicate completion
echo "Filtering complete. The file hipstr_ref_small.bed has been created."
