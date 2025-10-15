#!/bin/bash
# Script that finds the cancer variants that are specific to EUR and EAS populations and runs bedtools closest to find the unique variants that are closest to each other
cd ../preliminary_exploration/variant_selection/

# Input files
EUR_FILE="EUR_cancer_variants_genes.txt"
EAS_FILE="EAS_cancer_variants_genes.txt"

# Temporary files
EUR_SORTED="EUR_sorted.bed"
EAS_SORTED="EAS_sorted.bed"
EUR_SPECIFIC="EUR_specific.bed"
EAS_SPECIFIC="EAS_specific.bed"

# Output files
CLOSEST_OUTPUT="EUR_closest_EAS.txt"

# Step 0. Sort input files
sort -k1,1V -k2,2n $EUR_FILE > $EUR_SORTED
sort -k1,1V -k2,2n $EAS_FILE > $EAS_SORTED

# Step 1. Filter to population-specific variants

# Extract variant IDs from each file
cut -f4 $EUR_SORTED | sort | uniq > EUR_variant_ids.txt
cut -f4 $EAS_SORTED | sort | uniq > EAS_variant_ids.txt

#Find the variants that are unique to each population using the comm command
## the comm command compares sorted files FILE1 and FILE2 line by line. -2 suppress column 2 (lines unique to FILE2) -3 suppress column 3 (lines that appear in both files)
# Find EUR-specific variant IDs
comm -23 EUR_variant_ids.txt EAS_variant_ids.txt > EUR_specific_ids.txt
# Find EAS-specific variant IDs
comm -13 EUR_variant_ids.txt EAS_variant_ids.txt > EAS_specific_ids.txt

# Filter BED files to keep only specific variants
grep -F -f EUR_specific_ids.txt $EUR_SORTED > $EUR_SPECIFIC
grep -F -f EAS_specific_ids.txt $EAS_SORTED > $EAS_SPECIFIC

# Cleanup intermediate variant ID files
rm EUR_variant_ids.txt EAS_variant_ids.txt EUR_specific_ids.txt EAS_specific_ids.txt

# Step 2. Run bedtools closest
bedtools closest -a $EUR_SPECIFIC -b $EAS_SPECIFIC -d > temp_closest.txt

# Extract desired columns:
# EUR variant (cols 1-6), EAS variant (cols 7-12), distance (col 13)
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$8,$9,$10,$11,$12,$16,$17}' temp_closest.txt > $CLOSEST_OUTPUT

# Step 3. Cleanup
rm $EUR_SORTED $EAS_SORTED $EUR_SPECIFIC $EAS_SPECIFIC temp_closest.txt

echo "Population-specific closest variant analysis complete. Output: $CLOSEST_OUTPUT"

#Run bedtools intersect to find the genes that overlap with all cancer variants
bedtools intersect -a /exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/cancer_gwas.bed -b gencode_hg38_v47.gtf > all_cancer_genes.bed
