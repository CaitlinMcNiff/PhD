#!/bin/bash
## counts the number of eQTLs in each population file and outputs to a summary file

# Output file
output="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/GTEx_eQTLs/eQTL_counts.txt"
#> "$output"  # Empty the file if it already exists

# Loop through all population files (edit pattern if needed)
for file in /exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/GTEx_eQTLs/*_all_GTEx.bed; do
    # Extract population name (e.g., AFR from AFR_file.txt)
    population=$(basename "$file" | cut -d'_' -f1)

    # Count lines where the second-last column contains 'prom'
    total_count=$(wc -l < "$file")

    # Write result
    echo -e "${population}\t${total_count}" >> "$output"
done