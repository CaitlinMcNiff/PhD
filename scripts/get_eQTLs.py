import pandas as pd
import os
import glob
import subprocess

input_folder = '/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/GTEx_hg38_v10'  # Replace with your path
output_file = '/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/all_GTEx_hg38_v10.bed'
unique_lines = set()

# Find all matching parquet files
parquet_files = glob.glob(os.path.join(input_folder, "*.v10.eQTLs.signif_pairs.parquet"))

for parquet_file in parquet_files:
    print(f"Processing {os.path.basename(parquet_file)}")
    
    try:
        df = pd.read_parquet(parquet_file, columns=["variant_id"])
    except Exception as e:
        print(f"Error reading {parquet_file}: {e}")
        continue

    # Drop any missing or malformed variant IDs
    df = df.dropna(subset=["variant_id"])
    df = df[df["variant_id"].str.count("_") >= 4]  # ensure at least 4 underscores

    # Split variant_id into components
    parts = df["variant_id"].str.split("_", expand=True)
    df["chrom"] = parts[0]
    df["start"] = parts[1].astype(int)
    df["end"] = df["start"] + 1
    df["ref"] = parts[2]
    df["alt"] = parts[3]
    df["change"] = df["ref"] + "->" + df["alt"]
    df["dot"] = "."

    # Combine into tab-separated format
    for row in df.itertuples(index=False):
        line = f"{row.chrom}\t{row.start}\t{row.end}\t{row.variant_id}\t{row.change}\t{row.dot}"
        unique_lines.add(line)

# Write to file
with open(output_file, 'w') as f:
    for line in sorted(unique_lines):
        f.write(line + "\n")

print(f"Done. Output saved to: {output_file}")



# Run bedtool sintersect to find the number of phenotype-associated variants in each population that are associated with eQTL variants
subprocess.call("bedtools intersect -a gwas_1000_genomes/AFR_all_gwas.bed -b all_GTEx_hg38_v10.bed > AFR_all_GTEx.bed", shell=True)
subprocess.call("bedtools intersect -a gwas_1000_genomes/AMR_all_gwas.bed -b all_GTEx_hg38_v10.bed > AMR_all_GTEx.bed", shell=True)
subprocess.call("bedtools intersect -a gwas_1000_genomes/EAS_all_gwas.bed -b all_GTEx_hg38_v10.bed > EAS_all_GTEx.bed", shell=True)
subprocess.call("bedtools intersect -a gwas_1000_genomes/EUR_all_gwas.bed -b all_GTEx_hg38_v10.bed > EUR_all_GTEx.bed", shell=True)
subprocess.call("bedtools intersect -a gwas_1000_genomes/SAS_all_gwas.bed -b all_GTEx_hg38_v10.bed > SAS_all_GTEx.bed", shell=True)
print("bedtool sintersect run for each population with GETx file")