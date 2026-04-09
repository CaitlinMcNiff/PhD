#!/bin/bash

## Script that checks how many variants have an AF of 0 in all populations, and prints the number in variant_numbers.txt file

KGP_dir="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/1KGP_hg38" 
variant_nos="$KGP_dir/variant_numbers.txt" # output file

# Takes the 8th column (INFO) from the input file – which contains the AF data for all the populations, excluding any rows that start with #
#then splits columns by ; to get each population AF
awk -F'\t' '!/^#/ {print $8}' ERZ822766_merged.vcf | \
awk -F';' -v out="$variant_nos" '
BEGIN {
    count = 0 # initialise count
}
{
    delete info

    # loops over INFO fields - e.g. kv[1] = "EAS_AF", kv[2] = "0" -> info["EAS_AF"] = "0", info["EUR_AF"] = "0.2" etc.
    for (i = 1; i <= NF; i++) {
        split($i, kv, "=")
        info[kv[1]] = kv[2]
    }

    # checks AF is 0 for all 5 populations - treats 0, 0.0 & 0.00 all as zeros in case of different formatting
    if ((info["EAS_AF"] + 0) == 0 && 
        (info["EUR_AF"] + 0) == 0 &&
        (info["AFR_AF"] + 0) == 0 &&
        (info["AMR_AF"] + 0) == 0 &&
        (info["SAS_AF"] + 0) == 0) {
        count++
    }
}
END {
    print count
    print "Number of variants with AF=0 in all superpopulations -> " count >> out
}
'