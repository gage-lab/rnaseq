library(tidyverse)
library(tximport)
library(AnnotationDbi)
library(fishpond)

#' add feature ids to SummarizedExperiment-like object
#' @param obj SummarizedExperiment-like object
#' @param keys character vector of keys to map to feature IDs
#' @param gtf.df data.frame of gene info from GTF
#' @param gtf.cols character vector of columns from gtf.df to add to rowData
add_ids <- function(obj, keys, gtf.df, gtf.cols, orgdb) {
  # add gene names
  ids <- c("GENETYPE", "GENENAME", "SYMBOL", "ALIAS", "REFSEQ", "ENSEMBL", "ENTREZID")
  for (i in ids) {
    print(glue("Adding {i} feature IDs from org.Hs.eg.db"))
    rowData(obj)[[i]] <- AnnotationDbi::mapIds(orgdb, keys = keys, keytype = "ENSEMBL", column = i, multiVals = "first")
  }

  for (col in gtf.cols) {
    print(glue("Adding {col} from gtf"))
    rowData(obj)[[col]] <- gtf.df[rownames(obj), ][[col]]
  }
  return(obj)
}

#' Remove lowly expressed features
#' #' @param obj SummarizedExperiment-like object
#' @param minCount minimum number of counts for a feature to be kept
#' @param minN minimum number of samples with minCount for a feature to be kept
rm_lowexp <- function(obj, minCount = 10, minN = 2) {
  if (sum(grepl("^EN", rownames(obj))) < 10000) {
    # for testing
    minCount <- 1
  }
  # TODO: add condition variable for filtering
  print("Removing lowly expressed features")
  obj <- fishpond::labelKeep(obj, minCount = minCount, minN = minN)
  print(glue("Removing {sum(!rowData(obj)$keep)} features"))
  obj <- obj[mcols(obj)$keep, ]
  print(glue("{nrow(obj)} features remain"))
  return(obj)
}
