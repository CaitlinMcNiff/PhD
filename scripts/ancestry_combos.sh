#!/bin/bash

## This script counts the number of times a combination of populations was used in a study and how many studies used how many populations

DIR='/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS'

## Count the number of populations used in how many studies
outfile=${DIR}/study_population_counts.txt
echo "number_of_populations_studied   count" > ${outfile} # initialise file
awk -F'\t' '$14 != "" {print $3, $14}' ${DIR}/GWAS_ancestry_with_unique_populations.tsv \
| sort | uniq \
| awk '{print $1}' \
| sort | uniq -c \
| awk 'BEGIN {OFS="\t"}{counts[$1]++} END {for (c in counts) print c, counts[c]}' \
| sort -n >> ${outfile}

# Count the number of occurrences of different combinations of populations 
awk -F'\t' '
NR > 1 && $14 != "" {
    study = $3
    pop   = $14

    # makes a combined key for the study-population kair
    key = study SUBSEP pop
    if (!(key in seen)) {
        seen[key] = 1 # means that pair has already been seen
        pops[study] = pops[study] ? pops[study] "," pop : pop # stores the populations with that study
    }
}
END { # runs after code has gone through all file lines
    for (study in pops) {
        n = split(pops[study], arr, ",")# splits population string into array so populations can be manipulated individually

        # sort populations alphabetically within each study, so different populations orders are treated the same
        for (i = 1; i <= n; i++) {
            for (j = i + 1; j <= n; j++) {
                if (arr[i] > arr[j]) {
                    tmp = arr[i]
                    arr[i] = arr[j]
                    arr[j] = tmp
                }
            }
        }

        combo = arr[1]
        for (i = 2; i <= n; i++) combo = combo "," arr[i]

        count[combo]++
    }

    for (combo in count) {
        print count[combo], combo
    }
}' ${DIR}/GWAS_ancestry_with_unique_populations.tsv | sort -k1,1nr > ${DIR}/ancestry_combo_counts.txt