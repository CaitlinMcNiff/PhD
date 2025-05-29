setwd("/Volumes/phd")

library(dplyr)
library(readr)
library(stringr)

# Read ancestry file and rename column 2 as PMID
ancestry_df <- read_tsv("NHGRI_EBI_GWAS/GWAS_ancestry.tsv", col_names = TRUE)
colnames(ancestry_df)[3] <- "PMID"

# Define superpopulations and corresponding file names
superpops <- c("AFR", "AMR", "EUR", "EAS", "SAS")
file_paths <- paste0("gwas_1000_genomes/", superpops, "_all_gwas.bed")
variant_colnames <- c("CHR", "start_pos", "end_pos", "RSID", "Phenotype", "risk_AF", "p-value", "OR", "PMID", "v-CHR", "v-start_pos", "v-end_pos", "v-rsid", "REF", "ALT", "AF", "overlap")

# Initialize results table
results <- data.frame()

# Loop over each file/superpopulation
for (i in seq_along(superpops)) {
  pop <- superpops[i]
  file <- file_paths[i]
  
  # Read variant file
  variant_df <- read_tsv(file, col_names = FALSE)
  colnames(variant_df) <- variant_colnames
  
  # Extract PMIDs from ancestry that include this superpopulation
  studied_pmids <- ancestry_df %>%
    filter(grepl(pop, `SUPERPOPULATIONS`, ignore.case = TRUE)) %>%
    pull(PMID) %>%
    unique()
  
  # Count total, studied, and non-studied
  n_total <- nrow(variant_df)
  n_overlap <- sum(variant_df$PMID %in% studied_pmids)
  n_nonoverlap <- n_total - n_overlap
  
  # Append result
  results <- bind_rows(results, data.frame(
    Superpopulation = pop,
    Variants_Studied = n_overlap,
    Variants_Not_Studied = n_nonoverlap,
    Total_Variants = n_total
  ))
}

# Save summary
write_tsv(results, "gwas_ancestry_phenotype_overlap.tsv")
print(results)

results$p_value <- apply(results, 1, function(row) {
  a <- as.numeric(row["Variants_Studied"])
  b <- as.numeric(row["Variants_Not_Studied"])
  # Basic 2x2 test — assume same row/col structure
  mat <- matrix(c(a, b, sum(a, b) - a, sum(a, b) - b), nrow = 2)
  fisher.test(mat)$p.value
})
