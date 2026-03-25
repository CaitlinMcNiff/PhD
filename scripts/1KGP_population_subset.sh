#!/bin/bash
## Script to find the number of variants in the 1KGP data file and find the unmber of unique loci – allowing us to calculate how many loci appear more than once

KGP_dir="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/1KGP_hg38"

zcat $KGP_dir/1KGP_hg38/ERZ822766_merged.vcf.gz | sed '/^#/d' | awk '{FS="\t"} {print $1,$2,$4,$5}' > $KGP_dir/trimmed_variants.vcf

#count how many variants there are in the 1KGP data file
variant_nos="$KGP_dir/variant_numbers.txt"
total_vars=$(wc -l < $KGP_dir/trimmed_variants.vcf)
echo "total number of variants -> $total_vars" > "$variant_nos"

# sort file
sort -k1,1V -k2,2n $KGP_dir/trimmed_variants.vcf > $KGP_dir/1KGP_variants_sorted.bed 
echo "made sorted variant file"
# remove trimmed file – same as new file but not sorted.
rm $KGP_dir/trimmed_variants.vcf

# which loci appear more than once and how many are there
awk '{fs="\t"} {print $1,$2}' $KGP_dir/1KGP_variants_sorted.bed | uniq -c | sort -nr > $KGP_dir/variant_count.bed
echo "made a file with the loci variant counts"
uniq_vars=$(wc -l < $KGP_dir/1KGP_variants_sorted.bed)
echo "number of unique variants -> $uniq_vars" >> "$variant_nos"
echo "calculated unique variants"

# count how many unique loci there are – alternative way
#uniq_vars=$(awk -F'\t' '{print $1,$2}' $KGP_dir/1KGP_variants_sorted.bed | uniq |  wc -l)

# count how many loci have more than one variant
repeated_loci=$(awk -F'\t' '{print $1, $2}' "$KGP_dir/1KGP_variants_sorted.bed" | uniq -c | awk '$1 > 1' | wc -l)
echo "number of loci with multiple variants -> $repeated_loci" >> "$variant_nos"

# sanity check to make sure numbers are sensible
echo "TOTAL: $total_vars"
echo "UNIQUE: $uniq_vars"
echo "MULTIPLE: $repeated_loci"