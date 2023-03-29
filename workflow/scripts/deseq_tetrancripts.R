#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Jan 10, 2023

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(tximeta)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(glue)
})
snakemake@source("utilities.R")
options(readr.show_col_types = FALSE)


# make coldata for tximeta
coldata <- readr::read_tsv(snakemake@input$samplesheet)

files <- snakemake@input$quant
names(files) <- basename(dirname(dirname(files)))

coldata$files <- files[match(coldata$sample_name, names(files))]
coldata$names <- coldata$sample_name
coldata <- dplyr::select(coldata, -starts_with("fq"), -sample_name, -strandedness)

# remove low-quality samples
print(paste0("Removing samples: ", paste(snakemake@params[["samples_exclude"]], collapse = ", ")))
coldata <- dplyr::filter(coldata, !names %in% snakemake@params[["samples_exclude"]])

# get info from GTF
gene_gtf <- rtracklayer::readGFF(snakemake@input[["txome_gtf"]]) %>%
  dplyr::select(gene_id, gene_type, gene_name) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("gene_id")
rmsk_gtf <- rtracklayer::readGFF(snakemake@input[["rmsk_gtf"]]) %>%
  dplyr::select(gene_id, family_id, class_id) %>%
  dplyr::mutate(subfamily_id = gene_id, gene_id = paste(gene_id, family_id, class_id, sep = ":")) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("gene_id")

# DESeq
print("Reading TEtranscripts quantifications for DESeq")
se <- tximeta::tximeta(
  coldata,
  type = "none",
  txOut = TRUE,
  txIdCol = "gene_id",
  countsCol = "count",
  lengthCol = "length",
  abundanceCol = "TPM",
  importer = readr::read_tsv,
  skipMeta = TRUE
)

# add alternative names/ids
keys <- stringr::str_remove(rownames(rowData(se)), "\\..*$")
se <- add_ids(se, keys, gene_gtf, c("gene_name"))
se <- add_ids(se, keys, rmsk_gtf, c("class_id", "family_id", "subfamily_id"))
rowData(se)$gene_name[is.na(rowData(se)$gene_name)] <- rownames(rowData(se))[is.na(rowData(se)$gene_name)]

# Fit DESeq model
se <- rm_lowexp(se) # remove lowly expressed genes
form <- as.formula(snakemake@params$model)
dds <- DESeq2::DESeqDataSet(se, form)
print("Fitting DESeq model for differential gene expression (DGE) with transposable elements (TEs)")
dds <- DESeq2::DESeq(dds) # Fit DESeq model
saveRDS(dds, file = snakemake@output[[1]]) # save
