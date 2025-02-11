# Read and transform sample count data

library(tidyverse)
library(org.Hs.eg.db)

sample_count_file <- file.path("extdata", "GSE262138_Raw_count_data.txt")

dataset <- read_tsv(sample_count_file) |>
  mutate(Geneid = str_remove(Geneid, "\\.\\d+"))

mapped_gene_symbols <- select(org.Hs.eg.db, keys = dataset$Geneid, keytype = "ENSEMBL", columns = c("GENENAME", "SYMBOL"))

annotated_dataset <- dataset |>
  inner_join(mapped_gene_symbols, by = c(Geneid = "ENSEMBL")) |>
  filter(!is.na(SYMBOL)) |>
  summarise(across(where(is.numeric), ~ sum(.x)), .by = "SYMBOL") |>
  dplyr::select(Gene_Symbol = SYMBOL, starts_with("R")) |>
  write_csv(file.path("data", "GSE262138_Sample_Counts.csv"))
