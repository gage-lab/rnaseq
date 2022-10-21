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
read_tetranscripts <- function(counts, coldata, genes, rmsk, subfamily = "L1HS", deseq.obj = TRUE, 
		biotypes = c("protein_coding", "lncRNA", "antisense", "IG_C_gene", 
	"IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "IG_V_pseudogene", 
	"IG_J_pseudogene", "IG_C_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene",
	"TR_V_gene", "TR_V_pseudogene", "TR_J_pseudogene")){
    # biotypes selected from 10x genomics website
	# TODO: perform biotype filtering on the GTF before alignment

	# check inputs are strings
	for (input in list(counts, coldata, genes, rmsk)){
		stopifnot(class(input) == "character")
	}

	message(glue("Reading gene annotations from {genes}..."))
	message(glue("Using average length across basic transcripts for each gene..."))
    genes = rtracklayer::readGFF(genes) %>%
		dplyr::filter(tag == "basic", type == "transcript", transcript_biotype %in% biotypes) %>% # filter for basic transcripts of specified biotypes
		dplyr::mutate(length = end - start) %>%
		dplyr::group_by(gene_id, gene_name, transcript_biotype) %>% # group by genes for summarizing
		dplyr::summarize(length = mean(length), sd_length = sd(length)) %>% # take average length across transcripts for each gene
		dplyr::rename(feature_id = gene_id) %>%
		dplyr::select(feature_id, gene_name, length, transcript_biotype) 

	message(glue("Reading TE annotations from {rmsk}..."))
	rmsk = rtracklayer::readGFF(rmsk) %>%
		dplyr::mutate(feature_id = paste(gene_id, family_id, class_id, sep = ":"), biotype = "te") %>%
		dplyr::rename(subfamily_id = gene_id) %>%
		dplyr::select(feature_id, subfamily_id, family_id, class_id, biotype) %>%
		distinct()

	message(glue("Only including {subfamily} subfamily from TE annotations..."))
	if (subfamily == "L1HS"){
		rmsk = rmsk %>% 
			dplyr::filter(subfamily_id == subfamily) %>% 
			dplyr::mutate(length = 6000)
	}

	message(glue("Merging gene and TE annotations..."))
	annot = dplyr::bind_rows(genes, rmsk) %>% 
		dplyr::arrange(feature_id) %>%
		tibble::column_to_rownames("feature_id")

	message(glue("Reading sample metadata from {coldata}..."))
	samples = readr::read_tsv(coldata) %>% 
		column_to_rownames("sample_name") %>%
		dplyr::mutate_if(is.character, as.factor) 
    
	message(glue("Reading counts from {counts}..."))
    counts = readr::read_tsv(counts) %>% 
		dplyr::mutate(across(where(is.numeric), as.integer)) %>% # make sure all counts are integers  
		dplyr::arrange(feature_id) # sort by feature id

	message(glue("Removing {sum(!counts$feature_id %in% rownames(annot))}/{nrow(counts)} counted features not in annotation..."))
	counts = counts %>%
		dplyr::filter(feature_id %in% rownames(annot)) %>% # only include features in annotation
		tibble::column_to_rownames("feature_id") %>%
		as.matrix()

	stopifnot(dim(counts)[1] == dim(annot)[1])
	stopifnot(dim(counts)[2] == dim(samples)[1])

	# make matrix of gene lengths for each sample
	lengths = data.frame(row.names = rownames(annot))
	for (sample in colnames(counts)){
		lengths[sample] <- annot$length
	}

	message(glue("Creating DESeq2 object..."))
    # make SummarizedExperiment obj
    se = SummarizedExperiment(assays = SimpleList(counts = counts, lengths = lengths), colData = samples, rowData = annot)
	
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