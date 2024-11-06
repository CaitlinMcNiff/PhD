#!/bin/bash

# Define directories  
output_directory="gwas_1000_genomes" 
summary_file="$output_directory/bedtools_summary.txt"  # File to store line counts

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Initialize the summary file
echo "Intersect File Line Counts" > "$summary_file"
echo "==========================" >> "$summary_file"

# Step 1: Remove "chr" from the first column of all *_variants.vcf files
for vcf_file in *_variants.vcf; do
    # Remove "chr" from the first column 
    sed -i 's/^chr//' "$vcf_file" > "$output_directory/$(basename "$vcf_file" .vcf)_no_chr.vcf"
done

# Step 2: Run bedtools intersect for each combination of *_no_chr.vcf and *_gwas_*.bed files
for no_chr_vcf in "$output_directory"/*_no_chr.vcf; do
    for gwas_bed in *_gwas_*.bed; do
        # Extract base names for the output file
        vcf_base=$(basename "$no_chr_vcf" .vcf)
        bed_base=$(basename "$gwas_bed" .bed)
        
        # Output file name
        intersect_file="$output_directory/${vcf_base}_${bed_base}_intersect.bed"
        
        # Run bedtools intersect and output to a new file
        bedtools intersect -a "$no_chr_vcf" -b "$gwas_bed" > "$intersect_file"
        
        # Count the number of lines in the intersect file and add it to the summary
        line_count=$(wc -l < "$intersect_file")
        echo "${vcf_base}_vs_${bed_base}: $line_count lines" >> "$summary_file"
    done
done
