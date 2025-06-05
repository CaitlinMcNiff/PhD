import pandas as pd

# Load the data
file_path = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/gwas_catalog-ancestry_r2025-02-18.tsv"  
df = pd.read_csv(file_path, delimiter="\t")

# Define mapping of broad ancestry categories to 1KGP superpopulations
superpop_mapping = {
    "African American or Afro-Caribbean": "AFR",
    "African": "AFR",
    "African unspecified": "AFR",
    "Sub-Saharan African" : "AFR",
    "European": "EUR",
    "South Asian": "SAS",
    "East Asian": "EAS",
    "South East Asian": "EAS",  # Includes Malaysia, Singapore, Thailand, Vietnam, Indonesia (decision made using chatGPT)
    "Hispanic or Latin American": "AMR",
    "Admixed American": "AMR",
    "Native American": "AMR",
    "Asian unspecified": None,  # it is to broad to try to categorise it into either EAS or SAS with any confidence
    "NR": None  # NR (Not Reported) should be ignored
} # there are other broad ancestrys listed but they are too ambiguous or do not fit into a superpopulation, e.g. Middle Eastern or Oceanian

# Function to process ancestry column
def get_superpopulations(ancestry):
    #Splits multiple ancestries and maps them to 1000 Genomes superpopulations.
    if pd.isna(ancestry):
        return []
    ancestries = [a.strip() for a in ancestry.split(",")]  # Split by comma & remove spaces
    superpops = {superpop_mapping.get(a) for a in ancestries if a in superpop_mapping}  # Use a set to remove duplicates
    return [s for s in superpops if s]  # Remove None values

# Expands the DataFrame so that if a row contains a list of superpopulations, each item in the list will have its own row. This step is necessary to treat each superpopulation individually when counting them later
df["SUPERPOPULATIONS"] = df["BROAD ANCESTRAL CATEGORY"].apply(get_superpopulations)
df_exploded = df.explode("SUPERPOPULATIONS")  # Creates a new row for each superpopulation

# Drop duplicate entries for the same study + superpopulation
df_unique = df_exploded.drop_duplicates(subset=["STUDY ACCESSION", "SUPERPOPULATIONS"])

# Save the updated ancestry DataFrame to a new file
updated_file = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/GWAS_ancestry.tsv"
df_unique.to_csv(updated_file, sep="\t", header=True)
# Count unique studies per superpopulation
study_counts = df_unique.groupby("SUPERPOPULATIONS")["STUDY ACCESSION"].count()

# Save the results of the count data to a text file in the NHGRI_EBI_GWAS folder
output_file = "/exports/cmvm/datastore/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/NHGRI-EBI_ancestry.txt"
study_counts.to_csv(output_file, sep="\t", header=True)

print(f"Output saved to {output_file}")