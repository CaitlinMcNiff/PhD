#!/bin/bash

chromosome=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X')

output_dir="/exports/cmvm/datastore/sbms/groups/young-lab/caitlin/phd/1KGP_hg38"

for chr in "${chromosome[@]}"; do
    file="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    wget -P "$output_dir" "$file"
done