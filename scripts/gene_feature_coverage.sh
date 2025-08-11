#!/bin/bash

module load roslin/bedtools/2.31.1  # Load bedtools module

# === CONFIGURATION ===
GTF_FILE="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode/gencode_hg38_v47.gtf"  
GENOME_SIZE=3137300923  # GRCh38.p14 total non-N bases (from https://www.ncbi.nlm.nih.gov/grc/human/data)
features=("exon" "CDS" "transcript" "UTR" "start_codon" "stop_codon" "gene")  # Feature types to analyse
GEN_OUTFILE="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode/general_feature_coverage.tsv"

cd /exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode/ #|| exit  # Change to the directory containing the GTF file

# === HEADER ===
echo -e "Feature\tTotal_bp\tPercent_of_genome" > "${OUTFILE}"

# === MAIN LOOP ===
for FEATURE in "${features[@]}"; do
    echo "Processing ${FEATURE}..."

    # Extract lines matching the feature type
    awk -v feat="${FEATURE}" '$5 == feat' "${GTF_FILE}" > "${FEATURE}.gtf"
    # Sort and merge overlapping intervals – makes sure there are not ranges of a feature that overlap as this could inflate the total bp covered
    bedtools sort -i "${FEATURE}.gtf" | bedtools merge -i - > "${FEATURE}_merged.bed"
    # Calculate total bp covered
    TOTAL_BP=$(awk '{sum += $3 - $2} END {print sum}' "${FEATURE}_merged.bed")
    # Calculate % of genome
    PERCENT=$(echo "scale=6; ${TOTAL_BP} / ${GENOME_SIZE} * 100" | bc)

    # Write to output file
    echo -e "${FEATURE}\t${TOTAL_BP}\t${PERCENT}" >> "${GEN_OUTFILE}"
    # Clean up intermediate files
    rm "${FEATURE}.gtf" "${FEATURE}_merged.bed"
    done
echo "Coverage calculation complete. Results saved to ${GEN_OUTFILE}."


## CREATE THE SAME BUT FOR EACH POPULATION
# === CONFIGURATION ===
INPUT_FOLDER="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode/" 
POP_OUTPUT="population_feature_coverage1.tsv"

# === HEADER ===
echo -e "Population\tFeature\tTotal_bp\tPercent_of_genome" > "$POP_OUTPUT"

for file in "$INPUT_FOLDER"/*_gencode.vcf; do
    #filename=$(basename "$file")
    #population="${filename%%.*}"  # Extract population from filename, e.g., EUR from EUR.tsv
    population=$(basename "$file" "_gencode.vcf")
    
    echo "Processing population: $population"

    for feature in "${features[@]}"; do
        echo "  -> Feature: $feature"

        # Extract matching rows and build BED file from columns 1 (CHR), 2 (start), 3 (end)
        awk -v feat="$feature" 'BEGIN{FS="\t"; OFS="\t"} $10 == feat {print "chr"$1, $2, $3}' "$file" > "tmp_${population}_${feature}.bed"

        # Skip if file is empty
        if [[ ! -s tmp_${population}_${feature}.bed ]]; then
            echo "No data for feature: $feature"
            continue
        fi

        # Sort and merge overlapping intervals
        #bedtools sort -i "tmp_${population}_${feature}.bed" | bedtools merge -i - > "merged_${population}_${feature}.bed"

        # Calculate total base pairs
        #total_bp=$(awk '{sum += $3 - $2} END {print sum}' "merged_${population}_${feature}.bed")
        total_bp=$(sort "tmp_${population}_${feature}.bed" | uniq | wc -l)
        # Calculate percent coverage
        percent=$(echo "scale=6; $total_bp / $GENOME_SIZE * 100" | bc)

        # Write result
        echo -e "$population\t$feature\t$total_bp\t$percent" >> "$POP_OUTPUT"

        # Cleanup temp files
        #rm "tmp_${population}_${feature}.bed" "merged_${population}_${feature}.bed"
    done
done

echo "✅ Done! Results saved to: $POP_OUTPUT"