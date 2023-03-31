#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Dec 2, 2022

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(tximeta)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(DESeq2)
  library(fishpond)
  library(rtracklayer)
  library(glue)
})
snakemake@source("utilities.R")
options(readr.show_col_types = FALSE)

if (snakemake@config$species == "human") {
  orgdb <- org.Hs.eg.db
} else if (snakemake@config$species == "mouse") {
  orgdb <- org.Mm.eg.db
}

# make coldata for tximeta
coldata <- readr::read_tsv(snakemake@input$samplesheet)

files <- snakemake@input$quant
names(files) <- basename(dirname(dirname(files)))

coldata$files <- files[match(coldata$sample_name, names(files))]
coldata$names <- coldata$sample_name
coldata <- dplyr::select(coldata, -starts_with("fq"), -sample_name, -strandedness)

# remove low-quality samples
if (length(snakemake@config[["samples_exclude"]]) > 0) {
  stopifnot(all(snakemake@config[["samples_exclude"]] %in% coldata$names))
  print(paste0("Removing samples: ", paste(snakemake@config[["samples_exclude"]], collapse = ", ")))
  coldata <- dplyr::filter(coldata, !names %in% snakemake@config[["samples_exclude"]])
}

# get transcript and gene info from GTF
tx_gtf <- rtracklayer::readGFF(snakemake@input[["txome_gtf"]]) %>%
  dplyr::filter(!is.na(transcript_id)) %>%
  dplyr::select(gene_id, gene_type, gene_name, transcript_id, transcript_type, transcript_name, transcript_support_level) %>%
  dplyr::distinct()
gene_gtf <- dplyr::select(tx_gtf, gene_id, gene_type, gene_name) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("gene_id")
tx2gene <- dplyr::select(tx_gtf, gene_id, transcript_id) %>%
  dplyr::rename(GENEID = gene_id, TXNAME = transcript_id) %>%
  dplyr::relocate(TXNAME, GENEID)
rownames(tx_gtf) <- tx_gtf$transcript_id

# DESeq
print("Reading salmon quantifications for DESeq")
se <- tximeta::tximeta(coldata, txOut = FALSE, tx2gene = tx2gene, skipMeta = TRUE)
form <- as.formula(snakemake@params$model)
dds <- DESeq2::DESeqDataSet(se, form)
dds <- rm_lowexp(dds) # remove lowly expressed genes

# add alternative gene names/ids
keys <- stringr::str_remove(rownames(dds), "\\..*$")
dds <- add_ids(dds, keys, gene_gtf, c("gene_type", "gene_name"), orgdb = orgdb)

# fit DGE model
print("Fitting DESeq model for differential gene expression (DGE)")
dds <- DESeq2::DESeq(dds) # Fit DESeq model
saveRDS(dds, file = snakemake@output[["dge"]]) # save

# remove objects from memory
rm(se)
rm(dds)

# Swish
print("Reading salmon quantifications for Swish DTE/DTU")
se <- tximeta::tximeta(coldata, skipMeta = TRUE)
se <- fishpond::scaleInfReps(se)

# add alternative transcript names/ids
rowData(se)$gene_id <- tx2gene$GENEID[match(tx2gene$TXNAME, rownames(se))]
keys <- stringr::str_remove(rowData(se)$gene_id, "\\..*$")
se <- add_ids(se, keys, tx_gtf, c("gene_name", "transcript_type", "transcript_name", "transcript_support_level"), orgdb = orgdb)

se <- rm_lowexp(se) # remove lowly expressed transcripts

saveRDS(se, file = snakemake@output[["dte"]]) # save DTE
iso <- fishpond::isoformProportions(se)
saveRDS(iso, file = snakemake@output[["dtu"]]) # save DTU
