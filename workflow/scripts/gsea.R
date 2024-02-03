#!/usr/bin/env Rscript
# Author: Mike Cuoco, Lukas Karbacher
# Created on: Dec 12, 2022

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(zip)
  library(msigdbr)
  library(glue)
})
ggplot2::theme_set(ggplot2::theme_bw())
options(readr.show_col_types = FALSE)

# get inputs
dge <- readr::read_csv(snakemake@input[["dge"]])
gs <- snakemake@wildcards[["gs"]]

# # for debugging
# save.image(glue("{snakemake@wildcards[['gs']]}_{snakemake@wildcards[['contrast']]}_gsea.RData"))
# stop()

# make ranked list of genes
dge <- dplyr::filter(dge, !is.na(log2FoldChange))
ranked <- dge$log2FoldChange
names(ranked) <- dge$gene_id
ranked <- sort(ranked, decreasing = TRUE)

# get gene sets
gs_df <- readr::read_tsv(snakemake@input[["gs_df"]])
gs_list <- purrr::map(unique(gs_df$gs_name), function(gs) {
  dplyr::filter(gs_df, gs_name == gs) %>%
    dplyr::pull(human_ensembl_gene)
})
names(gs_list) <- unique(gs_df$gs_name)

# run GSEA
message("Running fgsea...")
res <- fgsea::fgsea(
  pathways = gs_list,
  stats = ranked
)

# save tabular results
readr::write_tsv(res, snakemake@output[[1]])
