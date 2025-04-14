library(dplyr)

# Define total length
total_length <- 78122275

# Define superpopulation counts
superpop_counts <- list(
  EAS = 9972280,
  AMR = 13494275,
  AFR = 19718874,
  EUR = 11155012,
  SAS = 12151736
)

# Function to perform bootstrapping and calculate 95% CI
bootstrap_mean <- function(count, total_length, n_bootstrap = 1000) {
  # Create binary vector
  gwas_vector <- c(rep(1, count), rep(0, total_length - count))
  
  # Perform bootstrapping and store means
  bootstrap_means <- replicate(n_bootstrap, mean(sample(gwas_vector, total_length, replace = TRUE)))
  
  # Calculate 95% confidence interval
  ci <- quantile(bootstrap_means, c(0.025, 0.975))
  
  return(c(mean(bootstrap_means), ci[1], ci[2])) # Return mean, lower, and upper CI
}

# Perform bootstrapping for each superpopulation and store results in a data frame
results <- do.call(rbind, lapply(names(superpop_counts), function(pop) {
  cat("Bootstrapping for", pop, "\n")
  c(pop, bootstrap_mean(superpop_counts[[pop]], total_length))
}))

# Convert results to a data frame
colnames(results) <- c("Superpopulation", "Mean", "Lower_CI", "Upper_CI")
results <- as.data.frame(results)

converted_results <- data.frame(Superpopulation = results$Superpopulation,
                                Lower_CI = as.numeric(results$Lower_CI) * total_length,
                                Upper_CI = as.numeric(results$Upper_CI) * total_length, 
                                Mean = as.numeric(results$Mean) * total_length)

# Save results to a text file
write.table(results, "raw_bootstrap_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(converted_results, "bootstrap_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Print results
print(results)

unique_counts <- list(
  EAS	= 1589651,
  AMR =	491541,
  AFR =	8378255,
  EUR	= 751881,
  SAS	= 1903891
)

# Perform bootstrapping for each unique variant superpopulation and store results in a data frame
results_unique <- do.call(rbind, lapply(names(superpop_counts), function(pop) {
  cat("Bootstrapping for", pop, "\n")
  c(pop, bootstrap_mean(unique_counts[[pop]], superpop_counts[[pop]]))
}))

colnames(results_unique) <- c("Superpopulation", "Mean", "Lower_CI", "Upper_CI")
results_unique <- as.data.frame(results_unique)

super_counts <- t(data.frame(superpop_counts))
super_counts <- as.data.frame(super_counts)
super_counts$Superpopulation <- rownames(super_counts)


converted_results_unique <- inner_join(results_unique, super_counts, by = "Superpopulation")

converted_results_unique <- converted_results_unique %>% mutate_at(c('Mean', 'Lower_CI', 'Upper_CI',  'V1'), as.numeric) %>% mutate(Mean = Mean * V1,
                                                                                                                                    Lower_CI = Lower_CI * V1, 
                                                                                                                                    Upper_CI = Upper_CI * V1)
converted_results_unique <- subset(converted_results_unique, select = -c(V1))

# Save results to a text file
write.table(results_unique, "raw_unique_variant_bootstrap_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(converted_results_unique, "unique_bootstrap_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)




##  Bootstrapping for the population-gwas overlap
# Line counts by phenotype group (overall total for each phenotype group)
phenotype_totals <- list(
  Cancer = 13580,
  Immunological = 41052,
  Neurological = 21626
)

# Input subgroup line counts
input_data <- data.frame(
  Superpopulation = c("AFR", "AFR", "AFR", "AMR", "AMR", "AMR", "EAS", "EAS", "EAS", 
                      "EUR", "EUR", "EUR", "SAS", "SAS", "SAS"),
  Phenotype = rep(c("Cancer", "Immunological", "Neurological"), times = 5),
  Line_Count = c(9458, 29182, 17834, 10479, 32386, 19164, 8909, 27482, 16862,
                 10593, 32564, 19275, 10009, 30996, 18617)
)

# Run bootstrap for each subgroup and store results
results <- do.call(rbind, lapply(1:nrow(input_data), function(i) {
   pop <- input_data$Superpopulation[i]
   phenotype <- input_data$Phenotype[i]
   count <- input_data$Line_Count[i]
   total <- phenotype_totals[[phenotype]]
  
   cat("Bootstrapping for", pop, "in", phenotype, "\n")
     
   boot_vals <- bootstrap_mean(count, total)
   data.frame(Superpopulation = pop,
              Phenotype = phenotype,
              Count = count,
              Mean = boot_vals[1],
              Lower_CI = boot_vals[2],
              Upper_CI = boot_vals[3])
 }))

pheno <- t(data.frame(phenotype_totals))
pheno <- as.data.frame(pheno)
pheno$Phenotype <- rownames(pheno)

converted_results_pheno <- inner_join(results, pheno, by = "Phenotype")
converted_results_pheno <- converted_results_pheno %>% mutate_at(c('Mean', 'Lower_CI', 'Upper_CI',  'V1'), as.numeric) %>% mutate(Mean = Mean * V1,
                                                                                                                                    Lower_CI = Lower_CI * V1, 
                                                                                                                                    Upper_CI = Upper_CI * V1)
converted_results_pheno <- subset(converted_results_pheno, select = -c(V1))

# Save results
write.table(results, file = "raw_bootstrap_superpop_phenotype_results_2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(converted_results_pheno, file = "bootstrap_superpop_phenotype_results_2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

print("Bootstrapping completed and results saved.")
