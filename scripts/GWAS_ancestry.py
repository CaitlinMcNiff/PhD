## This script maps the BROAD ANCESTRAL CATEGORY from the NHGRI-EBI GWAS Catalogue to 1KGP superpopulations, and counts the number of unique studies per superpopulation and stage.

import pandas as pd

# Load the data
file_path = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/gwas_catalog-ancestry_r2025-02-18.tsv"
df = pd.read_csv(file_path, delimiter="\t")

# Define mapping of broad ancestry categories to 1KGP superpopulations
superpop_mapping = {
    "African American or Afro-Caribbean": "AFR",
    "African": "AFR",
    "African unspecified": "AFR",
    "Sub-Saharan African": "AFR",
    "European": "EUR",
    "South Asian": "SAS",
    "East Asian": "EAS",
    "South East Asian": "EAS",  
    "Hispanic or Latin American": "AMR",
    "Admixed American": "AMR",
    "Native American": "AMR",
    "Asian unspecified": None,  
    "NR": None  
}

# Function to process ancestry column
def get_superpopulations(ancestry):
    if pd.isna(ancestry):
        return []
    ancestries = [a.strip() for a in ancestry.split(",")]  # Split by comma & remove spaces
    superpops = {superpop_mapping.get(a) for a in ancestries if a in superpop_mapping}  # Use a set to remove duplicates
    return [s for s in superpops if s]  # Remove None values

# Map to superpopulations
df["SUPERPOPULATIONS"] = df["BROAD ANCESTRAL CATEGORY"].apply(get_superpopulations)
df_exploded = df.explode("SUPERPOPULATIONS")  # Creates a new row for each superpopulation

# Drop duplicate entries for the same study + superpopulation
df_unique = df_exploded.drop_duplicates(subset=["STUDY ACCESSION", "SUPERPOPULATIONS", "STAGE"])

# Save the updated ancestry DataFrame
updated_file = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/GWAS_ancestry.tsv"
df_unique.to_csv(updated_file, sep="\t", header=True, index=False)

# ---- Counts ----
# Total studies per superpopulation
total_counts = df_unique.groupby("SUPERPOPULATIONS")["STUDY ACCESSION"].nunique()

# Counts per stage (Initial, Replication) for each superpopulation
stage_counts = (
    df_unique.groupby(["SUPERPOPULATIONS", "STAGE"])["STUDY ACCESSION"]
    .nunique()
    .unstack(fill_value=0)
)

# Merge total + stage counts
final_counts = pd.concat([total_counts, stage_counts], axis=1)
final_counts = final_counts.rename(columns={"STUDY ACCESSION": "Total_Studies"})

# Save results
output_file = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/NHGRI-EBI_ancestry.txt"
final_counts.to_csv(output_file, sep="\t", header=True)

print(f"Output saved to {output_file}")
