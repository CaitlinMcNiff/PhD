import pandas as pd
import re

# Define input and output file paths
input_file = "../NHGRI_EBI_GWAS/gwas_catalog_v1.0-associations_e113_r2025-02-18.tsv"
output_files = {
    "neurological": "../NHGRI_EBI_GWAS/neurological_gwas.bed",
    "immunological": "../NHGRI_EBI_GWAS/immunological_gwas.bed",
    "cancer": "../NHGRI_EBI_GWAS/cancer_gwas.bed"
}

# Define phenotype filters
filters = {
    "neurological": r"neuro|depression|brain|intel|schizo|bipolar|autis",
    "immunological": r"lymph|arthrit|inflam|immun|asthm|COVID|Treg|T\scell|B\scell|white\sblood|neutrophil|basophil",
    "cancer": r"cancer"
}

# Define exclusion terms for cancer
cancer_exclusion = r"excluded|Non-c|Illnesses"

# Load the file
df = pd.read_csv(input_file, delimiter="\t", dtype=str)  # Load all as strings

# Extract and clean chromosome information
def get_chromosome(chrom_col, pos_col):
    """Returns chromosome from column 12 unless missing, then extracts from column 22."""
    if pd.isna(chrom_col) or chrom_col == "NA":
        match = re.match(r'chr(\d+|X|Y|MT):\d+', str(pos_col))
        return match.group(1) if match else None
    if "," in chrom_col:  # Ignore rows with multiple chromosome values
        return None
    return chrom_col

# Process the data
df["chromosome"] = df.apply(lambda row: get_chromosome(row.iloc[11], row.iloc[21]), axis=1)
df["position"] = pd.to_numeric(df.iloc[:, 12], errors="coerce")  # Convert position column to numeric (NaN if invalid)
df["position_plus1"] = df["position"].fillna(0).astype(int) + 1  # Fill NaNs with 0, convert to int, then add 1
df["position"] =df["position"].astype('Int64') # converts the position column to an integer – Int64 accounts for NA values
df["phenotype"] = df.iloc[:, 7]
df["risk_AF"] = df.iloc[:, 26]
df["p_value"] = df.iloc[:, 27]
df["OR"] = df.iloc[:, 30]
df["rsid"] = df.iloc[:, 21]  # Only include rsID if chromosome exists

# Keep only rows with a valid chromosome
df = df.dropna(subset=["chromosome"])

# Select relevant columns
df = df[["chromosome", "position", "position_plus1", "rsid", "phenotype", "risk_AF", "p_value", "OR"]]

# Function to filter and save data
def filter_and_save(df, category, regex, exclusion=None):
    filtered_df = df[df["phenotype"].str.contains(regex, case=False, na=False, regex=True)]
    if exclusion:
        filtered_df = filtered_df[~filtered_df["phenotype"].str.contains(exclusion, case=False, na=False, regex=True)]
    
    # Sort by chromosome (natural order) and position (numeric order)
    filtered_df = filtered_df.sort_values(by=["chromosome", "position"], key=lambda x: pd.to_numeric(x, errors="coerce"))

    # Save to file without header/index
    filtered_df.to_csv(output_files[category], sep="\t", index=False, header=False)

# Generate output files
for category, regex in filters.items():
    exclude = cancer_exclusion if category == "cancer" else None
    filter_and_save(df, category, regex, exclude)

print("Filtered GWAS files have been saved successfully!")
