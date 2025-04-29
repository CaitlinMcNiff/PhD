import gzip
import os

output_file = gzip.open('/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/gencode_hg38_v47.gtf.gz', 'wb')

with gzip.open("/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/gencode.v47.annotation.gtf.gz", 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue

        cols = line.strip().split('\t')
        chr = cols[0].replace('chr', '')
        feature = cols[2]
        start_pos = cols[3]
        end_pos = cols[4]
        strand = cols[6]
        info = cols[8]

        output_file.write(f"{chr}\t{start_pos}\t{end_pos}\t{strand}\t{feature}\t{info}\n")

output_file.close()






