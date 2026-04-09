#!/bin/bash

KGP_dir="/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/1KGP_hg38"
variant_nos="$KGP_dir/variant_numbers.txt"

awk -F'\t' '!/^#/ {print $8}' ERZ822766_merged.vcf | \
awk -F';' -v out="$variant_nos" '
BEGIN {
    count = 0
}
{
    delete info

    for (i = 1; i <= NF; i++) {
        split($i, kv, "=")
        info[kv[1]] = kv[2]
    }

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