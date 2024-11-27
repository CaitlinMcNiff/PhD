#!/bin/bash

# load bedtools module
module load igmm/apps/BEDTools/2.31.1

output_directory="gwas_1000_genomes"
summary_file="$output_directory/bedtools_summary.txt"

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Initialize the summary file
echo "Intersect File Line Counts" > "$summary_file"
echo "==========================" >> "$summary_file"

# Create the files for the subgroups of the openGWAS file  – immunological, neurological, and cancer
#grep 'neuro\|depression\|brain\|intel\|schizo\|bipolar\|autism' /exports/cmvm/datastore/sbms/groups/young-lab/human/openGWAS_genome_wide_sig.bed > neurological_gwas.bed
#grep 'lymph\|arthrit\|inflam\|immun\|asthm\|COVID\|Treg\|T\scell\|B\scell\|white\sblood\|neutrophil\|basophil' /exports/cmvm/datastore/sbms/groups/young-lab/human/openGWAS_genome_wide_sig.bed > immunological_gwas.bed 
#grep 'cancer' /exports/cmvm/datastore/sbms/groups/young-lab/human/openGWAS_genome_wide_sig.bed | grep -v 'excluded\|Non-c\|Illnesses' > cancer_gwas.bed

# Remove "chr" from the first column of all the openGWAS files so they can be matched to superpopulation files
for gwas_file in *_gwas.bed; do
    sed -i 's/^chr//' "$gwas_file"
done

# Run bedtools intersect for each combination of *_variants.vcf and *_gwas_*.bed files
for pop_vcf in *_variants_sorted.bed; do
    for gwas_file in *_gwas.bed; do
        # extract base names for the output file
        vcf_base=$(basename "$pop_vcf" _variants_sorted.bed)
        bed_base=$(basename "$gwas_file" .bed)

        #output file name
        intersect_file=$output_directory/${vcf_base}_${bed_base}.bed

        bedtools intersect -a "$pop_vcf" -b "$gwas_file" > "$intersect_file"

        # Count the number of lines in the intersect file and add it to the summary
        line_count=$(wc -l < "$intersect_file")
        echo "${vcf_base}_vs_${bed_base}: $line_count lines" >> "$summary_file"
    done
    # bedtools intersect for each of the superpopulations for all GWAS hits
    oG_file=$output_directory/${vcf_base}_openGWAS.bed
    bedtools intersect -a /exports/cmvm/datastore/sbms/groups/young-lab/human/openGWAS_genome_wide_sig.bed -b "$pop_vcf" > "$oG_file"

done