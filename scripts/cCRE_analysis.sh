#!/bin/bash
## counts the number of cCREs in each population file and outputs to a summary file

# Output file
output="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/cCRE/ccre_counts.txt"
> "$output"  # Empty the file if it already exists
echo -e "Population\tpromoter\tenhancer\ttotal_ccres" > "${output}"

# Loop through all population files (edit pattern if needed)
for file in /exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/cCRE/*_all_ccres.bed; do
    # Extract population name (e.g., AFR from AFR_file.txt)
    population=$(basename "$file" | cut -d'_' -f1)

    # Count lines where the second-last column contains 'prom'
    count_prom=$(awk -F'\t' '{ if ($(NF-1) ~ /prom/) ++c } END { print c+0 }' "$file")
    count_enh=$(awk -F'\t' '{ if ($(NF-1) ~ /enh/) ++c } END { print c+0 }' "$file")
    total_count=$(wc -l < "$file")

    # Write result
    echo -e "${population}\t${count_prom}\t${count_enh}\t${total_count}" >> "$output"
done