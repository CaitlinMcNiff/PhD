import gzip
import os
import subprocess

# --- 1. Process and simplify the GENCODE GTF file ---
gencode_input_path = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/gencode.v47.annotation.gtf.gz"
gencode_output_path = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary/gencode_hg38_v47.gtf.gz"

with gzip.open(gencode_input_path, 'rt') as infile, gzip.open(gencode_output_path, 'wb') as outfile:
    for line in infile:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        if len(cols) < 9:
            continue  # skip malformed lines
        chr_clean = cols[0].replace('chr', '')
        feature = cols[2]
        start_pos = cols[3]
        end_pos = cols[4]
        strand = cols[6]
        info = cols[8]

        output_line = f"{chr_clean}\t{start_pos}\t{end_pos}\t{strand}\t{feature}\t{info}\n"
        outfile.write(output_line.encode('utf-8'))

# --- 2. Process population-specific GWAS files into simplified BED files ---

populations = ['EAS', 'AMR', 'AFR', 'EUR', 'SAS']

# Input/output directories
gwas_base_path = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd"
bed_dir = "gwas_1000_genomes"
output_dir = os.path.join(gwas_base_path, "preliminary")

# Column names in GWAS files (no headers)
variant_colnames = [
    "CHR", "start_pos", "end_pos", "RSID", "Phenotype",
    "risk_AF", "p_value", "OR", "PMID", "v_CHR",
    "v_start_pos", "v_end_pos", "v_rsid", "REF", "ALT", "AF", "overlap"
]

for pop in populations:
    input_path = os.path.join(bed_dir, f"{pop}_all_gwas.bed")
    output_path = os.path.join(output_dir, f"{pop}_all_gwas.bed")

    with open(input_path, 'rt') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            cols = line.strip().split('\t')
            if len(cols) < 8:
                continue  # skip malformed lines

            chr = cols[0]
            start_pos = cols[1]
            end_pos = cols[2]
            rsid = cols[3]
            pheno = cols[4]

            outfile.write(f"{chr}\t{start_pos}\t{end_pos}\t{rsid}\t{pheno}\n")

# --- 3. Run bedtools intersect for each population ---

gtf_unzipped = gencode_output_path.rstrip('.gz')  # Expected output filename after unzip
gtf_for_bedtools = gtf_unzipped

# Unzip the GTF if needed (bedtools needs plain text input)
if not os.path.exists(gtf_unzipped):
    subprocess.run(["gunzip", "-k", gencode_output_path])

#header = [
    "CHR", "start_pos", "end_pos", "RSID", "Phenotype",          # from BED
    "chr_gtf", "gtf_start", "gtf_end", "gtf_strand", "gtf_feature", "gtf_info",  # from GTF
    "overlap_length"  # from -wo]

# Run bedtools intersect
for pop in populations:
    input_bed = os.path.join(output_dir, f"{pop}_all_gwas.bed")
    output_vcf = os.path.join(output_dir, f"{pop}_gencode.vcf")

    command = [
        "bedtools", "intersect",
        "-a", input_bed,
        "-b", gtf_for_bedtools,
        "-wo"
    ]

    with open(output_vcf, "w") as outfile:
        subprocess.run(command, stdout=outfile)

    print(f"Finished bedtools for {pop}: {output_vcf}")