#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Feb 6, 2023

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})
options(readr.show_col_types = FALSE)

res <- readr::read_csv(snakemake@input[[1]])

# for debugging
# save.image(glue("{snakemake@wildcards[['contrast']]}_ora.RData"))

# UP-regulated
sig <- res %>%
  dplyr::filter(log2FoldChange > snakemake@config$de$cutoffs$log2FoldChange, padj < snakemake@config$de$cutoffs$FDR) %>%
  dplyr::pull("gene_id") %>%
  stringr::str_remove(".\\d+$")
sigFC <- res %>%
  dplyr::filter(log2FoldChange > snakemake@config$de$cutoffs$log2FoldChange, padj < snakemake@config$de$cutoffs$FDR) %>%
  dplyr::pull("log2FoldChange")

names(sigFC) <- sig

# run GO enrichment
message(glue("Running GO enrichment for {length(sig)} DOWN-regulated genes..."))
ego <- clusterProfiler::enrichGO(
  gene = sig,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  universe = stringr::str_remove(res$gene_id, ".\\d+$"),
  readable = TRUE,
  keyType = "ENSEMBL",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
)

# save results
as.data.frame(ego) %>% readr::write_tsv(snakemake@output$up)

# DOWN-regulated
sig <- res %>%
  dplyr::filter(log2FoldChange < -snakemake@config$de$cutoffs$log2FoldChange, padj < snakemake@config$de$cutoffs$FDR) %>%
  dplyr::pull("gene_id") %>%
  stringr::str_remove(".\\d+$")
sigFC <- res %>%
  dplyr::filter(log2FoldChange < -snakemake@config$de$cutoffs$log2FoldChange, padj < snakemake@config$de$cutoffs$FDR) %>%
  dplyr::pull("log2FoldChange")

# run GO enrichment
message(glue("Running GO enrichment for {length(sig)} UP-regulated genes..."))
ego <- clusterProfiler::enrichGO(
  gene = sig,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  universe = stringr::str_remove(res$gene_id, ".\\d+$"),
  readable = TRUE,
  keyType = "ENSEMBL",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
)

# save results
as.data.frame(ego) %>% readr::write_tsv(snakemake@output$down)
