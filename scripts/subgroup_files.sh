#!/bin/bash

# Create the files for the subgroups of the NHGRI-EBI GWAS catalogue, selecting for the columns that contain the chr, pos, rsid, phenotype, risk AF, p-val, and OR/BETA & sorted by chr and pos
# neurological phenotypes
awk 'BEGIN {FS="\t"}; {print $12, $13, $22, $8, $27, $28, $31}' ../NHGRI_EBI_GWAS/gwas_catalog_v1.0-associations_e113_r2025-02-18.tsv | grep 'neuro\|depression\|brain\|intel\|schizo\|bipolar\|autis' | sort -k1,1V -k2,2n > ../NHGRI_EBI_GWAS/neurological_gwas.bed
#immunological phenotypes
awk 'BEGIN {FS="\t"}; {print $12, $13, $22, $8, $27, $28, $31}' ../NHGRI_EBI_GWAS/gwas_catalog_v1.0-associations_e113_r2025-02-18.tsv | grep 'lymph\|arthrit\|inflam\|immun\|asthm\|COVID\|Treg\|T\scell\|B\scell\|white\sblood\|neutrophil\|basophil'| sort -k1,1V -k2,2n > ../NHGRI_EBI_GWAS/immunologicl_gwas.bed
# cancer phenotypes
awk 'BEGIN {FS="\t"}; {print $12, $13, $22, $8, $27, $28, $31}' ../NHGRI_EBI_GWAS/gwas_catalog_v1.0-associations_e113_r2025-02-18.tsv | grep 'cancer' /exports/cmvm/datastore/sbms/groups/young-lab/human/openGWAS_genome_wide_sig.bed | grep -v 'excluded\|Non-c\|Illnesses' | sort -k1,1V -k2,2n > ../NHGRI_EBI_GWAS/cancer_gwas.bed