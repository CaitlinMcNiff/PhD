#!/bin/bash

# Define input files
EUR_VARIANTS="../gwas_1000_genomes/EUR_cancer_gwas.bed"
EAS_VARIANTS="../gwas_1000_genomes/EAS_cancer_gwas.bed"
GENCODE_GENES="../preliminary_exploration/gencode/gencode_hg38_v47.gtf"

# Define output files
EUR_intersect_output="../preliminary_exploration/variant_selection/EUR_cancer_variants_genes.txt"
EAS_intersect_output="../preliminary_exploration/variant_selection/EAS_cancer_variants_genes.txt"

# Ensure files are sorted for bedtools
sort -k1,1V -k2,2n $EUR_VARIANTS > EUR_sorted.bed
sort -k1,1V -k2,2n $EAS_VARIANTS > EAS_sorted.bed
sort -k1,1V -k2,2n $GENCODE_GENES > gencode_sorted.bed

# Intersect EUR cancer variants with GENCODE genes
bedtools intersect -a EUR_sorted.bed -b gencode_sorted.bed -wa -wb > temp_EUR_genes.txt

# Select desired columns (chr, start, end, variant_id, gene_id)
awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"$4"\t"$23}' temp_EUR_genes.txt | awk 'BEGIN {FS = ";"} {print $1, $4}' | grep gene_name | uniq > $EUR_intersect_output

# Intersect EAS cancer variants with GENCODE genes
bedtools intersect -a EAS_sorted.bed -b gencode_sorted.bed -wa -wb > temp_EAS_genes.txt

# Select desired columns
awk 'BEGIN {FS="\t"}{print $1"\t"$2"\t"$3"\t"$4"\t"$23}' temp_EAS_genes.txt | awk 'BEGIN {FS = ";"} {print $1, $4}' | grep gene_name | uniq > $EAS_intersect_output

# Cleanup temporary files
rm temp_EUR_genes.txt temp_EAS_genes.txt EUR_sorted.bed EAS_sorted.bed gencode_sorted.bed

echo "Intersection complete. Outputs:"
echo "$EUR_intersect_output"
echo "$EAS_intersect_output"

# Tidy up the file bedtools file that has all cancer variants
awk 'BEGIN {OFS="\t"}{FS="\t"} {print $1, $2, $3, $4, $5, $15}' all_cancer_genes.bed | awk 'BEGIN {FS = ";"} {print $1, $4}' | grep gene_name > all_cancer_varants_genes.txt
# Count occurrences of each gene and sort by frequency
awk 'BEGIN {FS=" "}{OFS="\t"} {print $NF}' all_cancer_variants_genes.txt | sort | uniq -c | sort -k1,1r > cancer_gene_count.txt