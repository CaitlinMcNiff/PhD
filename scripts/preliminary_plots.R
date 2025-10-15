# Scripts to make the figures that show the difference between the different superpopulation variant numbers
library(dplyr)
library(ggplot2)
library(tidyverse)

setwd("/Volumes/cmvm/sbms/groups/young-lab/caitlin/phd/") # working directory when connected to the datatore folder on laptop

# read in files for the total number of varaints for each superpopulation
pop <- read.table("1KGP_hg38/pop_count_hg38.txt", sep = "\t")
colnames(pop) <- c("Superpopulation", "Total_Lines")
bootstrapped <- read.table("bootstrap_results.txt", sep = "\t", header = T) # bootstrapping results for the total variant count for
pop_boot <- full_join(pop, bootstrapped, by = "Superpopulation")
# adjusting the lines so they are not in scientific notation
pop_boot$Adjusted_lines <- pop_boot$Total_Lines / 1e7
pop_boot$Adjusted_uCI <- pop_boot$Upper_CI / 1e7
pop_boot$Adjusted_lCI <- pop_boot$Lower_CI / 1e7


# read in files for the  number of varaints unique to each superpopulation
unique_pop <- read.table("1KGP_hg38/unique_pop_count_hg38.txt", sep = "\t")
colnames(unique_pop) <- c("Superpopulation", "Lines")
unique_bootstrap <- read.table("unique_bootstrap_results.txt", sep = "\t", header = T)
uniq_pop_boot <- full_join(unique_pop, unique_bootstrap, by = "Superpopulation")
# adjusting the lines so they are not in scientific notation
uniq_pop_boot$Adjusted_lines <- uniq_pop_boot$Lines / 1e6
uniq_pop_boot$Adjusted_uCI <- uniq_pop_boot$Upper_CI / 1e6
uniq_pop_boot$Adjusted_lCI <- uniq_pop_boot$Lower_CI / 1e6

# plotting total variants per population
pop_plot <- ggplot(pop_boot, aes(Superpopulation, Adjusted_lines, fill = Superpopulation)) + geom_bar(stat="identity") +
  geom_errorbar(aes(x = Superpopulation, ymin = Adjusted_lCI, ymax = Adjusted_uCI)) +
  ylab(bquote("Number of Variants "(x10^7))) +
  xlab("Superpopulation") +
  theme_classic() +
  theme(legend.position="none")
pop_plot
# Save high-resolution image
ggsave("population_variant_count.png", pop_plot, dpi = 300, width = 6, height = 4)

# plotting the number of unique variant per population
unique_pop_plot <- ggplot(uniq_pop_boot, aes(Superpopulation, Adjusted_lines, fill = Superpopulation)) + geom_bar(stat="identity") +
  geom_errorbar(aes(x = Superpopulation, ymin = Adjusted_lCI, ymax = Adjusted_uCI)) +
  ylab(bquote("Number of Unique Variants "(x10^6))) +
  xlab("Superpopulation") +
  theme_classic() +
  theme(legend.position="none")
unique_pop_plot
ggsave("unique_population_variant_count.png", unique_pop_plot, dpi = 300, width = 6, height = 4)


## bedtools summary of the overlap between superpopulations and phenotype groups
# reading in the bedtools summary file and manipulating the rows so they can be plotted
lines <- read_lines("gwas_1000_genomes/bedtools_summary.txt")[-1] # reads in the lines of the file, excluding the first line
bedtools <- tibble(raw = lines) %>%
  mutate(raw = str_remove(raw, " lines$")) %>%  # Remove " lines" from the end
  separate(raw, into = c("Superpopulation", "Phenotype"), sep = "_vs_", extra = "merge") %>%  # Split at "_vs_"
  separate(Phenotype, into = c("Phenotype", "Line_Count"), sep = ": ") %>%  # Split at ": "
  mutate(Phenotype = str_remove(Phenotype, "_gwas$"),  # Remove "_gwas"
         Line_Count = as.integer(Line_Count))  # Convert count to integer

bed_boot <- read.table("bootstrap_superpop_phenotype_results.txt", sep = "\t", header = T)
bedtools_boot <- full_join(bedtools, bed_boot, by = c("Superpopulation", "Phenotype"))
bedtools_boot <- subset(bedtools_boot, select = -c(Line_Count))
# creating a new dataframe that merges the dataframes with the number of variants and the phenotypes so the proportion can be calculated
proportion_df <- full_join(pop, bedtools_boot, by ="Superpopulation")
proportion_df <- proportion_df[order(proportion_df$Superpopulation),]
proportion_df$Proportion <- proportion_df$Count / proportion_df$Total_Lines *100
proportion_df$Adjusted_uCI <- proportion_df$Upper_CI / proportion_df$Total_Lines *100
proportion_df$Adjusted_lCI <- proportion_df$Lower_CI / proportion_df$Total_Lines *100

phenotype_order <-  c('Cancer', 'Neurological', 'Immunological') # changing the order of the phenotypes so it is more aesthetic
# plotting the raw numbers of the overlap between the variants and the phenotype groups
bedtools_plot <- ggplot(bedtools_boot, aes(factor(Phenotype, level = phenotype_order), Count, fill = Superpopulation)) + 
  geom_bar(stat = "identity",position=position_dodge(), color="white",) + 
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), 
                position = position_dodge(width = 0.9), 
                width = 0.25, 
                color = "black") +
  xlab("Phenotype") +
  ylab("Number of Variants") +
  theme_classic()
bedtools_plot
ggsave("phenotype_variant_count.png", bedtools_plot, dpi = 300, width = 7, height = 5)


# plotting the proportion of the phenotype overlap compared to the total number of variants for the population
prop_plot <- ggplot(proportion_df, aes(factor(Phenotype, level = phenotype_order), Proportion, fill = Superpopulation)) + 
  geom_bar(stat = "identity",position=position_dodge(), color="white",) + 
  geom_errorbar(aes(ymin = Adjusted_lCI, ymax = Adjusted_uCI), 
                position = position_dodge(width = 0.9), 
                width = 0.25, 
                color = "black") +
  xlab("Phenotype") +
  ylab("Propotion of the total variants (%)") +
  theme_classic()
prop_plot
ggsave("phenotype_variant_proportion.png", prop_plot, dpi = 300, width = 7, height = 5)

