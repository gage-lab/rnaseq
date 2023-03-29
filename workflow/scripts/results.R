#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Dec 12, 2022

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(fishpond)
  library(glue)
  library(tidyverse)
})

se <- readRDS(snakemake@input[[1]])
coef <- snakemake@wildcards$contrast

terms <- as.formula(snakemake@params$model) %>%
  terms() %>%
  labels()

stopifnot(length(terms) >= 1)

condition <- stringr::str_split_1(coef, "_")[1]
condition1 <- stringr::str_split_1(coef, "_vs_")[1] %>%
  stringr::str_remove_all(paste0("^", condition, "_"))
condition2 <- stringr::str_split_1(coef, "_vs_")[2]

contrast <- c(condition, condition1, condition2)
if (snakemake@wildcards[["de"]] %in% c("dge", "dge_te")) {
  # get results
  print(glue("Getting DESeq2 results for {condition}: {condition1} vs {condition2}"))
  res <- DESeq2::results(se, contrast = contrast, cooksCutoff = FALSE)

  # perform LFC shrinkage
  print("Performing LFC shrinkage")
  res <- DESeq2::lfcShrink(dds = se, coef = coef, type = "apeglm", res = res)

  # add annotations
  df <- rowData(se) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::select(-tidyselect::contains("_vs_"))

  res <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(df)
} else if (snakemake@wildcards[["de"]] %in% c("dte", "dtu")) {
  # run swish for DTE/DTU
  print(glue("Running Swish for {condition}: {condition1} vs {condition2} for {toupper(snakemake@wildcards[['de']])} analysis"))
  se <- se[, se[[condition]] %in% c(condition1, condition2)] # only keep two conditions for test
  se$condition <- factor(se[[condition]], levels = c(condition1, condition2)) # set condition as factor

  # run swish
  if (length(terms) == 2) {
    dte <- fishpond::swish(y = se, x = terms[1], cov = terms[2])
  } else if (length(terms) == 1) {
    dte <- fishpond::swish(y = se, x = terms[1])
  } else {
    stop("Only one or two terms are allowed in the model formula when using swish.")
  }
  res <- rowData(dte) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("tx_id") %>%
    dplyr::mutate(padj = qvalue, log2FoldChange = log2FC, baseMean = 10^log10mean)
}

# save results
print(paste0("Saving results to ", snakemake@output[[1]]))
readr::write_csv(res, file = snakemake@output[[1]])
