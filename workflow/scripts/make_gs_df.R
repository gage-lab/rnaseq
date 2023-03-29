#!/usr/bin/env Rscript
# Author: Mike Cuoco, Lukas Karbacher, Thais Sabedot
# Created on: Dec 13, 2022

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
    library(readr)
    library(msigdbr)
})
options(readr.show_col_types = FALSE)

# save.image(file = glue("{snakemake@wildcards[["gs"]]}.RData"))

if (snakemake@wildcards[["gs"]] == "hallmark") {
    gs.df <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
} else if (snakemake@wildcards[["gs"]] == "kegg") {
    gs.df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
} else if (snakemake@wildcards[["gs"]] == "tft") {
    gs.df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
} else if (snakemake@wildcards[["gs"]] == "go") {
    gs.df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
    gs.df <- rbind(gs.df, msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC"))
    gs.df <- rbind(gs.df, msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF"))
}

readr::write_tsv(gs.df, snakemake@output[[1]])
