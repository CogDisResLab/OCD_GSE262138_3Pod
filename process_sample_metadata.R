# Extract the sample metadata from the Series Matrix File

library(tidyverse)

series_matrix_file <- "extdata/GSE262138_series_matrix.txt"

dataset <- read_tsv(series_matrix_file, skip = 31, col_names = FALSE) |>
    mutate(rowname = str_remove(X1, "!") |> str_c("x", row_number())) |>
    filter(
        str_detect(rowname, "data_processing", negate = TRUE),
        str_detect(rowname, "matrix_table", negate = TRUE)
    ) |>
    column_to_rownames() |>
    t() |>
    as_tibble() |>
    filter(str_detect(Sample_titlex1, "Sample", negate = TRUE)) |>
    mutate(
        Tissue = str_extract(Sample_characteristics_ch1x10, "tissue: (\\w+)", 1),
        CellType = str_extract(Sample_characteristics_ch1x11, "cell type: (.*)", 1),
        Sex = str_extract(Sample_characteristics_ch1x12, "Sex: ([M|F])", 1),
        Diagnosis = str_extract(Sample_characteristics_ch1x13, "group/disease: (.*)", 1),
        RIN = str_extract(Sample_characteristics_ch1x14, "rin: (.*)", 1) |> as.numeric() |> round(2),
        SampleID = str_extract(Sample_titlex1, "R\\d+")
    ) |>
    select(
        GEOID = Sample_geo_accessionx2,
        SampleID,
        Organism = Sample_organism_ch1x9,
        Tissue, CellType, Sex, RIN, Diagnosis
    ) |>
    mutate(Diagnosis = if_else(Diagnosis == "control", "CTRL", Diagnosis)) |>
    write_csv(file.path("data", "GSE262138_Sample_Metadata.csv"))
