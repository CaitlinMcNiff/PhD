import pandas as pd

# Load the data
file_path = "/exports/cmvm/datastore/sbms/groups/young-lab/caitlin/phd/NHGRI_EBI_GWAS/gwas_catalog-ancestry_r2025-02-18.tsv"  # Update with your file path
df = pd.read_csv(file_path, delimiter="\t")

# Define mapping of broad ancestry categories to 1KGP superpopulations
superpop_mapping = {
    "African American or Afro-Caribbean": "AFR",
    "African": "AFR",
    "European": "EUR",
    "South Asian": "SAS",
    "East Asian": "EAS",
    "Hispanic or Latin American": "AMR",
    "Admixed American": "AMR"
}

# Apply mapping to create a new column for superpopulation
df["SUPERPOPULATION"] = df["BROAD ANCESTRAL CATEGORY"].map(superpop_mapping)

# Count unique GWAS studies per superpopulation
study_counts = df.groupby("SUPERPOPULATION")["STUDY ACCESSION"].nunique()

# Display results
print(study_counts)
