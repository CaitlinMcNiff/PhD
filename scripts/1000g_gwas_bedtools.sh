#!/bin/bash

# load bedtools module
module load igmm/apps/BEDTools/2.31.1

KGP_dir="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/1KGP_hg38"
gwas_dir="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS"
output_directory="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/gwas_1000_genomes"
summary_file="$output_directory/bedtools_summary.txt"

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Initialize the summary file
echo "Intersect File Line Counts" > "$summary_file"

# Run bedtools intersect for each combination of *_variants.vcf and *_gwas_*.bed files
for pop_vcf in "$KGP_dir"/*_variants_sorted.vcf; do
    for gwas_file in "$gwas_dir"/*_gwas.bed; do
        # extract base names for the output file
        vcf_base=$(basename "$pop_vcf" _variants_sorted.vcf)
        bed_base=$(basename "$gwas_file" .bed)

        #output file name
        intersect_file=$output_directory/${vcf_base}_${bed_base}.bed

        bedtools intersect -a "${gwas_file}" -b "${pop_vcf}" -wo > "${intersect_file}"

        # Count the number of lines in the intersect file and add it to the summary
        line_count=$(wc -l < "$intersect_file")
        echo "${vcf_base}_vs_${bed_base}: ${line_count} lines" >> "${summary_file}"
    done
    
done