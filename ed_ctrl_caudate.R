# ED v CTRL DLPFC

library(tidyverse)
library(edgeR)
source("process_dge_fun.R")

counts_file <- file.path("data", "GSE262138_Sample_Counts.csv")
metadata_file <- file.path("data", "GSE262138_Sample_Metadata.csv")

metadata <- read_csv(metadata_file) |>
  filter(Diagnosis %in% c("CTRL", "ED"), Tissue != "DLPFC")

counts <- read_csv(counts_file) |>
  column_to_rownames("Gene_Symbol") |>
  select(any_of(metadata$SampleID))


groups <- metadata$Diagnosis
design_matrix <- model.matrix(~ 0 + groups)

dge <- DGEList(counts = counts, group = groups)

keep <- filterByExpr(dge)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

out_filtered <- process_dge(dge_filtered, design_matrix, num_genes = 1000)

out_filtered$complete_table |> 
  rownames_to_column("gene_name") |> 
  select(gene_name, logFC, P.Value = PValue) |>
  write_csv("results/ed_vs_control_caudate.csv")

save(out_filtered, file = "results/ed_vs_control_caudate.RData")
