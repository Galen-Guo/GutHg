library(maaslin3)
library(tidyverse)
# Input data table ####

### Read abundance table ####
# taxa_table_name <- system.file("extdata", "/UOttawa/Metagenomics/work/taxonomy_species.tsv", package = "maaslin3")
taxa_table <- read_tsv("/UOttawa/Metagenomics/work/taxonomy_species.tsv") %>%
  mutate(taxon = str_replace(taxon," ", "_")) %>%
  column_to_rownames("taxon") %>%
  t() %>% as.data.frame()

### Read metadata table ####
# metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package = "maaslin3")
metadata <- read_tsv("/UOttawa/Metagenomics/work/metadata.tsv") %>%
  mutate(Individual =  factor(Individual),
         Treatment  = factor(Treatment),
         Time = factor(Time)) %>%
  column_to_rownames("SampleID")

#Running MaAsLin 3 ####

set.seed(1)
fit_out <- maaslin3(input_data = as.data.frame(taxa_table),
                    input_metadata = as.data.frame(metadata),
                    output = 'hmp2_output',
                    formula = '~ Individual  + Time',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1,
                    plot_summary_plot = TRUE,
                    plot_associations = TRUE)
