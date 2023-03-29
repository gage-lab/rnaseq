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

res <- readr::read_csv(snakemake@input$dge)

#' Run GO analysis
#' @param res DESeq2 results
#' @param up_or_down "UP" or "DOWN"
#' @param results_fn filename for results
myGO <- function(res, up_or_down, results_fn) {
    # error checking
    stopifnot(up_or_down %in% c("UP", "DOWN"))

    if (up_or_down == "UP") {
        # UP-regulated
        sig <- res %>%
            dplyr::filter(log2FoldChange > snakemake@params[["LFCcutoff"]], padj < snakemake@params[["FDRcutoff"]]) %>%
            dplyr::pull("gene_id") %>%
            stringr::str_remove(".\\d+$")
        sigFC <- res %>%
            dplyr::filter(log2FoldChange > snakemake@params[["LFCcutoff"]], padj < snakemake@params[["FDRcutoff"]]) %>%
            dplyr::pull("log2FoldChange")
    } else {
        # DOWN-regulated
        sig <- res %>%
            dplyr::filter(log2FoldChange < -snakemake@params[["LFCcutoff"]], padj < snakemake@params[["FDRcutoff"]]) %>%
            dplyr::pull("gene_id") %>%
            stringr::str_remove(".\\d+$")
        sigFC <- res %>%
            dplyr::filter(log2FoldChange < -snakemake@params[["LFCcutoff"]], padj < snakemake@params[["FDRcutoff"]]) %>%
            dplyr::pull("log2FoldChange")
    }

    names(sigFC) <- sig

    # run GO enrichment
    message(glue("Running GO enrichment for {length(sig)} {up_or_down}-regulated genes..."))
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
    as.data.frame(ego) %>% readr::write_tsv(results_fn)
}

myGO(res, up_or_down = "UP", results_fn = snakemake@output$"resultsUP")
myGO(res, up_or_down = "DOWN", results_fn = snakemake@output$"resultsDOWN")
