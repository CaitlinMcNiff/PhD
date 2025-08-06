#!/bin/bash

# === CONFIGURATION ===
GTF_FILE="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode/gencode_hg38_v47.gtf"  
GENOME_SIZE=3088269832                  # GRCh38 total base pairs (excluding Ns)
FEATURES=("exon" "CDS" "transcript" "UTR" "start_codon" "stop_codon" "gene")  # Feature types to analyse
OUTFILE="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode/feature_coverage.tsv"

cd /exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/gencode/

# === HEADER ===
echo -e "Feature\tTotal_bp\tPercent_of_genome" > "$OUTFILE"

# === MAIN LOOP ===
for FEATURE in "${FEATURES[@]}"; do
    echo "Processing ${FEATURE}..."

    # Extract lines matching the feature type
    awk -v feat="${FEATURE}" '$5 == feat' "${GTF_FILE}" > "${FEATURE}.gtf"

    # Sort and merge overlapping intervals
    bedtools sort -i "${FEATURE}.gtf" | bedtools merge -i - > "${FEATURE}_merged.bed"

    # Calculate total bp covered
    TOTAL_BP=$(awk '{sum += $3 - $2} END {print sum}' "${FEATURE}_merged.bed")

    # Calculate % of genome
    PERCENT=$(echo "scale=6; ${TOTAL_BP} / ${GENOME_SIZE} * 100" | bc)

    # Write to output file
    echo -e "${FEATURE}\t${TOTAL_BP}\t${PERCENT}" >> "${OUTFILE}"

    # Clean up intermediate files (optional)
    rm "${FEATURE}.gtf" "${FEATURE}_merged.bed"
done

echo "Coverage calculation complete. Results saved to ${OUTFILE}."
