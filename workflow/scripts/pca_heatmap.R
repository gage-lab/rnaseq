#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Dec 6, 2022

suppressPackageStartupMessages({
    library(tidyverse)
    library(PCAtools)
    library(DESeq2)
    library(paletteer)
    library(ComplexHeatmap)
})

# get inputs/outputs
con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")
dds <- readRDS(snakemake@input[[1]])
form <- as.formula(snakemake@params$model)
terms <- form %>%
    terms() %>%
    labels()

# variance stabilizing transformation
print("Performing variance stabilizing transformation")
vst <- DESeq2::rlog(dds, blind = TRUE)

# PCA
# TODO: add loadings plot/output
print("Running PCA")
p <- PCAtools::pca(assay(vst), metadata = colData(dds))
bi <- PCAtools::biplot(p, legendPosition = "right", lab = NULL) +
    coord_fixed()
ggsave(plot = bi, filename = snakemake@output$pca_plot, width = 10, height = 10)
scree <- PCAtools::screeplot(p) + theme_classic()
ggsave(plot = scree, filename = snakemake@output$scree_plot, width = 10, height = 10)
comps <- getComponents(p, seq_len(5))
pairs <- PCAtools::pairsplot(p, components = comps[!is.na(comps)], gridlines.major = FALSE, gridlines.minor = FALSE) + coord_fixed()
ggsave(plot = pairs, filename = snakemake@output$pairs_plot)


# Heatmaps
plot_df <- colData(dds) %>%
    as.data.frame() %>%
    dplyr::select(all_of(c(labels(terms(form)))))
col_ha <- ComplexHeatmap::HeatmapAnnotation(df = plot_df)
row_ha <- ComplexHeatmap::HeatmapAnnotation(df = plot_df, which = "row", show_legend = FALSE, show_annotation_name = FALSE)

# euclidean distance
dst <- t(assay(vst)) %>%
    dist(diag = TRUE, upper = TRUE) %>%
    as.matrix()
svg(filename = snakemake@output$dist_heatmap, width = 10, height = 8)
p <- ComplexHeatmap::Heatmap(dst, bottom_annotation = col_ha, right_annotation = row_ha, name = "euclidean dst", column_title = "Euclidean Distance")
print(p)
dev.off()

# pearson correlation
corr <- cor(assay(vst))
svg(filename = snakemake@output$corr_heatmap, width = 10, height = 8)
p <- ComplexHeatmap::Heatmap(corr, bottom_annotation = col_ha, right_annotation = row_ha, name = "pearson", column_title = "Pearson Correlation")
print(p)
dev.off()
