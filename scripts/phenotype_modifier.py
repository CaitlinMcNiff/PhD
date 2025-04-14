import pandas as pd
import re

# Define input and output file paths
input_file = "../NHGRI_EBI_GWAS/gwas_catalog_v1.0-associations_e113_r2025-02-18.tsv"
output_files = {
    "neurological": "../NHGRI_EBI_GWAS/neurological_gwas.bed",
    "immunological": "../NHGRI_EBI_GWAS/immunological_gwas.bed",
    "cancer": "../NHGRI_EBI_GWAS/cancer_gwas.bed",
    "all": "../NHGRI_EBI_GWAS/all_gwas.bed"
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

# Function to extract chromosome and position
def extract_chrom_and_pos(chrom_col, pos_col):
    """Returns a tuple (chromosome, position). Extracts chromosome from column 12 unless missing, 
    then extracts it from column 21 and retrieves position from the format 'chrN:position'."""
    
    if pd.isna(chrom_col) or chrom_col == "NA":
        match = re.match(r'chr(\d+|X|Y|MT):(\d+)', str(pos_col))  # Extract chromosome and position
        if match:
            return match.group(1), int(match.group(2))  # Chromosome, Position
        return None, None  # If no match, return None
    
    if "," in chrom_col:  # Ignore rows with multiple chromosome values
        return None, None
    
    if ";" in chrom_col:  # Ignore rows with multiple chromosome values
        return None, None
    
    return chrom_col, None  # Return chromosome from col 12, but keep position unchanged

# Apply function to extract chromosome and position
df[["chromosome", "new_position"]] = df.apply(
    lambda row: pd.Series(extract_chrom_and_pos(row.iloc[11], row.iloc[21])),
    axis=1
)

# Process position column
df["position"] = pd.to_numeric(df.iloc[:, 12], errors="coerce")  # Convert column 13 to numeric
df["position"] = df["new_position"].combine_first(df["position"])  # Use extracted position if available
df["position_plus1"] = df["position"].fillna(0).astype(int) + 1  # Handle NaNs and add 1
df["position"] = df["position"].astype('Int64')  # Keep NA values

# Extract relevant columns
df["pmid"] = df.iloc[:,1]
df["phenotype"] = df.iloc[:, 7]
df["risk_AF"] = df.iloc[:, 26]
df["p_value"] = df.iloc[:, 27]
df["OR"] = df.iloc[:, 30]

# If chromosome was extracted from column 21, set rsid to empty string
df["rsid"] = df.apply(lambda row: "" if pd.isna(row.iloc[11]) or row.iloc[11] == "NA" else row.iloc[21], axis=1)

# Replace missing values in the final three columns with "NR"
df[["risk_AF", "p_value", "OR", "rsid", "pmid"]] = df[["risk_AF", "p_value", "OR", "rsid", "pmid"]].fillna("NR")

# Keep only rows with a valid chromosome
df = df.dropna(subset=["chromosome"])

# Select relevant columns
df = df[["chromosome", "position", "position_plus1", "rsid", "phenotype", "risk_AF", "p_value", "OR", "pmid"]]

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

df_sorted = df.sort_values(by=["chromosome", "position"], key=lambda x: pd.to_numeric(x, errors="coerce"))
df_sorted.to_csv(output_files["all"], sep="\t", index=False, header=False)

print("Filtered GWAS files have been saved successfully!")
