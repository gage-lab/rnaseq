#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Oct 20, 2022
#
# Create a DESeq2 object from TEcount output

# Load packages and Set Options ------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(glue))
options(readr.num_threads = snakemake@threads, readr.show_col_types = FALSE)

#' Create a DESeq2 object from
#' TEcount .cntTable file, sample metadata, TE annotation, and gene annotation
#' 
#' @param coldata, a data.frame with the following columns: files, names, ...
#' @param deseq.obj, a logical to indicate whether to return data in DESeq2 format
read_tetranscripts <- function(counts, coldata, genes, rmsk, deseq.obj = TRUE){
    
	# check inputs are strings
	for (input in list(counts, coldata, genes, rmsk)){
		stopifnot(class(input) == "character")
	}

	message(glue("Reading gene annotations from {genes}..."))
    genes = rtracklayer::readGFF(genes) %>%
		dplyr::filter(type == "gene") %>%
		dplyr::rename(id = gene_id, biotype = gene_biotype) %>%
		dplyr::select(id, gene_name, biotype) 

	message(glue("Reading TE annotations from {rmsk}..."))
	rmsk = rtracklayer::readGFF(rmsk) %>%
		dplyr::mutate(id = paste(gene_id, family_id, class_id, sep = ":"), biotype = "te") %>%
		dplyr::rename(subfamily_id = gene_id) %>%
		dplyr::select(id, subfamily_id, family_id, class_id, biotype) %>%
		distinct()
	annot = dplyr::bind_rows(genes, rmsk) %>% 
		dplyr::arrange(id) %>%
		tibble::column_to_rownames("id")

	message(glue("Reading sample metadata from {coldata}..."))
	samples = readr::read_tsv(coldata) %>% 
		column_to_rownames("sample_name") %>%
		dplyr::mutate_if(is.character, as.factor) 
    
	message(glue("Reading counts from {counts}..."))
	# make sure all counts are integers  
    counts = readr::read_tsv(counts) %>% 
		dplyr::mutate(across(where(is.numeric), as.integer)) %>%
		dplyr::arrange(id) %>%
		tibble::column_to_rownames("id") %>%
		as.matrix()

    # make RangedSummarizedExperiment obj
    se = SummarizedExperiment(assays = SimpleList(counts = counts), colData = samples, rowData = annot)
	
	# make DESeq2 obj
	dds = DESeqDataSet(se, design = ~ 1)
	return(dds)
}

dds = read_tetranscripts(counts = snakemake@input[["counts"]], 
					rmsk = snakemake@input[["rmsk"]], 
					genes = snakemake@input[["genes"]],
					coldata = snakemake@input[["coldata"]])

message(glue("Saving DESeq object to {snakemake@output[[1]]}..."))
saveRDS(dds, snakemake@output[[1]])