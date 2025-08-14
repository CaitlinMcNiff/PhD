#!/bin/bash
# Get the breakdown of feature types from the GENCODE GTF file
awk 'BEGIN {FS = "\t"}{print $5}' gencode_hg38_v47.gtf | sort | uniq -c | sort -nr > gencode_region_breakdown.txt

# Get the breakdown of gene types from the GENCODE GTF file – filters by only the lines that contain 'gene' as theif feature type
awk 'BEGIN {FS = "\t"}{print $5, $6}' gencode_hg38_v47.gtf | grep -w gene | awk 'BEGIN {FS = ";"} {print $2}' | sort | uniq -c | sort -nr  > gencode_gene_types.txt

output_dir="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode" 
summary_file1="$output_dir/gencode_population_region_summary.txt"  # File to store line counts
summary_file2="$output_dir/gencode_population_gene_summary.txt"

# Clear existing summary files
> "$summary_file1"
> "$summary_file2"

# Loop through all population-specific intersect files
for file in "$output_dir"/*_gencode.vcf; do
    pop_name=$(basename "$file" "_gencode.vcf")

    # Get the breakdown of feature types for the population
    echo "===== $pop_name =====" >> "$summary_file1"
    awk 'BEGIN {FS = "\t"}{print $10}' "$file" | sort | uniq -c | sort -nr >> "$summary_file1"
    echo "" >> "$summary_file1"

    # Get the breakdown of gene types from the population
    echo "===== $pop_name =====" >> "$summary_file2"
    awk 'BEGIN {FS = "\t"}{print $10, $11}' "$file" | grep -w gene | awk 'BEGIN {FS = ";"} {print $2}' | sort | uniq -c | sort -nr >> "$summary_file2"
    echo "" >> "$summary_file2"

    # Filter the output files from the gencode analysis so there is only one line for each variant and corresponding gene ID
    awk 'BEGIN {FS = "\t"}{print $1, $2, $3, $11}' "$file" | awk 'BEGIN {FS = ";"} {print $1}' | sort | uniq > ${pop_name}_gencode_uniq.bed
done