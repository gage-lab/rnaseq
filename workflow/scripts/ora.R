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
  library(org.Mm.eg.db)
})
options(readr.show_col_types = FALSE)

res <- readr::read_csv(snakemake@input[[1]])

if (snakemake@config$species == "human") {
  orgdb <- org.Hs.eg.db
} else if (snakemake@config$species == "mouse") {
  orgdb <- org.Mm.eg.db
}

#' Run GO analysis
#' @param res DESeq2 results
#' @param up_or_down "UP" or "DOWN"
#' @param results_fn filename for results
#' @param plot_fn filename for plot pdf
myGO <- function(res, up_or_down, results_fn, plot_fn) {
  # error checking
  stopifnot(up_or_down %in% c("UP", "DOWN"))

  if (up_or_down == "UP") {
    # UP-regulated
    sig <- res %>%
      dplyr::filter(log2FoldChange > snakemake@config[["de"]][["cutoffs"]][["log2FoldChange"]], padj < snakemake@config[["de"]][["cutoffs"]][["FDR"]]) %>%
      dplyr::pull("gene_id")
    sigFC <- res %>%
      dplyr::filter(log2FoldChange > snakemake@config[["de"]][["cutoffs"]][["log2FoldChange"]], padj < snakemake@config[["de"]][["cutoffs"]][["FDR"]]) %>%
      dplyr::pull("log2FoldChange")
  } else {
    # DOWN-regulated
    sig <- res %>%
      dplyr::filter(log2FoldChange < -snakemake@config[["de"]][["cutoffs"]][["log2FoldChange"]], padj < snakemake@config[["de"]][["cutoffs"]][["FDR"]]) %>%
      dplyr::pull("gene_id")
    sigFC <- res %>%
      dplyr::filter(log2FoldChange < -snakemake@config[["de"]][["cutoffs"]][["log2FoldChange"]], padj < snakemake@config[["de"]][["cutoffs"]][["FDR"]]) %>%
      dplyr::pull("log2FoldChange")
  }

  names(sigFC) <- sig

  # run GO enrichment
  message(glue("Running GO enrichment for {length(sig)} {up_or_down}-regulated genes..."))
  ego <- clusterProfiler::enrichGO(
    gene = sig,
    OrgDb = orgdb,
    ont = "ALL",
    universe = res$gene_id,
    readable = TRUE,
    keyType = "ENSEMBL",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
  )

  message(glue("Found {nrow(ego)} enriched GO terms for {length(sig)} {up_or_down}-regulated genes"))

  # save results
  as.data.frame(ego) %>% readr::write_tsv(results_fn)

  ego <- clusterProfiler::enrichGO(
    gene = sig,
    OrgDb = orgdb,
    ont = "ALL",
    universe = res$gene_id,
    readable = TRUE,
    keyType = "ENSEMBL",
    pvalueCutoff = 1,
    qvalueCutoff = 0.05,
  )

  if (length(ego) == 0 || nrow(ego) <= 3) {
    file.create(plot_fn)
    return(warning("Less than 4 enriched GO terms found, writing empty plots"))
  }

  # plot results
  ego2 <- enrichplot::pairwise_termsim(ego)

  pdf(plot_fn, width = 14, height = 10)
  p <- enrichplot::treeplot(ego2)
  print(p)
  p <- enrichplot::cnetplot(ego, categorySize = "pvalue", foldChange = sigFC)
  print(p)
  p <- enrichplot::emapplot(ego2)
  print(p)
  dev.off()
}

myGO(res, up_or_down = "UP", results_fn = snakemake@output[["resultsUP"]], plot_fn = snakemake@output[["plotUP"]])
myGO(res, up_or_down = "DOWN", results_fn = snakemake@output[["resultsDOWN"]], plot_fn = snakemake@output[["plotDOWN"]])
