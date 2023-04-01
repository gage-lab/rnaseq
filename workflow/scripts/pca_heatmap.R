#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Dec 6, 2022

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(PCAtools)
  library(DESeq2)
  library(paletteer)
  library(ComplexHeatmap)
})

# get inputs/outputs
dds <- readRDS(snakemake@input[[1]])
form <- as.formula(snakemake@params$model)
terms <- form %>%
  terms() %>%
  labels()

# variance stabilizing transformation
message("Performing regularized log transformation")
vst <- DESeq2::rlog(dds, blind = TRUE)

# PCA
# TODO: add loadings plot/output
message("Running PCA")
p <- PCAtools::pca(assay(vst), metadata = colData(dds))

# make plots
pdf(snakemake@output$pca)
scree <- PCAtools::screeplot(p) + theme_classic()
print(scree)

comps <- getComponents(p, seq_len(5))
for (t in terms) {
  bi <- PCAtools::biplot(p, legendPosition = "right", colby = t, lab = NULL) +
    coord_fixed()
  print(bi)
  pairs <- PCAtools::pairsplot(p, components = comps[!is.na(comps)], gridlines.major = FALSE, gridlines.minor = FALSE) + coord_fixed()
  print(pairs)
}
dev.off()

# Heatmaps
plot_df <- colData(dds) %>%
  as.data.frame() %>%
  dplyr::select(all_of(c(terms)))
col_ha <- ComplexHeatmap::HeatmapAnnotation(df = plot_df)
row_ha <- ComplexHeatmap::HeatmapAnnotation(df = plot_df, which = "row", show_legend = FALSE, show_annotation_name = FALSE)

# euclidean distance
pdf(snakemake@output$heatmap)
dst <- t(assay(vst)) %>%
  dist(diag = TRUE, upper = TRUE) %>%
  as.matrix()
p <- ComplexHeatmap::Heatmap(dst, bottom_annotation = col_ha, right_annotation = row_ha, name = "euclidean dst", column_title = "Euclidean Distance")
print(p)

# pearson correlation
corr <- cor(assay(vst))
p <- ComplexHeatmap::Heatmap(corr, bottom_annotation = col_ha, right_annotation = row_ha, name = "pearson", column_title = "Pearson Correlation")
print(p)
dev.off()
