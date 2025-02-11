# Read and transform sample count data

library(tidyverse)

sample_count_file <- file.path("extdata", "GSE262138_Raw_count_data.txt")

dataset <- read_tsv(sample_count_file) |>
    write_csv(file.path("data", "GSE262138_Sample_Counts.csv"))
