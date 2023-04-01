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
  library(org.Mm.eg.db)
  library(DESeq2)
  library(glue)
})
snakemake@source("utilities.R")
options(readr.show_col_types = FALSE)

# set orgdb
if (snakemake@config$species == "human") {
  orgdb <- org.Hs.eg.db
} else if (snakemake@config$species == "mouse") {
  orgdb <- org.Mm.eg.db
}

# handle output of TEtranscripts vs TElocal
if (snakemake@wildcards$quant_level == "subfamily") {
  rmsk_gtf <- rtracklayer::readGFF(snakemake@input[["rmsk_gtf"]]) %>%
    dplyr::select(gene_id, family_id, class_id) %>%
    dplyr::mutate(subfamily_id = gene_id, gene_id = paste(gene_id, family_id, class_id, sep = ":")) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames("gene_id")
} else if (snakemake@wildcards$quant_level == "locus") {
  rmsk_gtf <- rtracklayer::readGFF(snakemake@input[["rmsk_gtf"]]) %>%
    dplyr::select(transcript_id, gene_id, family_id, class_id) %>%
    dplyr::mutate(locus_id = transcript_id, subfamily_id = gene_id, gene_id = paste(transcript_id, gene_id, family_id, class_id, sep = ":")) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames("gene_id")
}

# make coldata for tximeta
coldata <- readr::read_tsv(snakemake@input$samplesheet)

files <- snakemake@input$quant
names(files) <- basename(dirname(files))

coldata$files <- files[match(coldata$sample_name, names(files))]
coldata$names <- coldata$sample_name
coldata <- dplyr::select(coldata, -starts_with("fq"), -sample_name)

# remove low-quality samples
if (length(snakemake@config[["samples_exclude"]]) > 0) {
  stopifnot(all(snakemake@config[["samples_exclude"]] %in% coldata$names))
  message(paste0("Removing samples: ", paste(snakemake@config[["samples_exclude"]], collapse = ", ")))
  coldata <- dplyr::filter(coldata, !names %in% snakemake@config[["samples_exclude"]])
}

# get info from GTF
gene_gtf <- rtracklayer::readGFF(snakemake@input[["txome_gtf"]]) %>%
  dplyr::select(gene_id, gene_type, gene_name) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("gene_id")

# DESeq
message("Reading TE quantifications for DESeq")
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

# remove lowly expressed genes/TEs
se <- rm_lowexp(se)

# add alternative names/ids
keys <- stringr::str_remove(rownames(rowData(se)), "\\..*$")
se <- add_ids(se, keys, gene_gtf, c("gene_name"), orgdb = orgdb)

rowData(se)[["class_id"]] <- rmsk_gtf[rownames(se), ][["class_id"]]
rowData(se)[["family_id"]] <- rmsk_gtf[rownames(se), ][["family_id"]]

if (snakemake@wildcards$quant_level == "locus") {
  rowData(se)[["locus_id"]] <- rmsk_gtf[rownames(se), ][["locus_id"]]
}

rowData(se)$gene_name[is.na(rowData(se)$gene_name)] <- rownames(rowData(se))[is.na(rowData(se)$gene_name)]

# Fit DESeq model
form <- as.formula(snakemake@params$model)
dds <- DESeq2::DESeqDataSet(se, form)
message("Fitting DESeq model for differential gene expression (DGE) with transposable elements (TEs)")
dds <- DESeq2::DESeq(dds) # Fit DESeq model
saveRDS(dds, file = snakemake@output[[1]]) # save
