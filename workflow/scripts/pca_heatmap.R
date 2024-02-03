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
  library(glue)
  library(circlize)
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
nsamples <- ncol(assay(vst))

# PCA
# TODO: add loadings plot/output
message("Running PCA")
p <- PCAtools::pca(assay(vst), metadata = colData(dds))

# make plots
pdf(snakemake@output$pca, width = 10, height = 10)
ncomp <- min(10, ncol(assay(vst)))
comps <- getComponents(p, seq_len(ncomp))
scree <- PCAtools::screeplot(p, components = comps) + theme_classic()
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

# # create color mappings
map_colors <- function(column_name) {
  # check if numeric
  if (snakemake@config$de$datatypes[[column_name]] == "categorical") {
    col <- paletteer::paletteer_d("ggsci::nrc_npg", n = length(unique(plot_df[[column_name]]))) %>% as.vector()
    names(col) <- levels(factor(plot_df[[column_name]]))
  } else if (snakemake@config$de$datatypes[[column_name]] == "numeric") {
    # make continuous color palette
    breaks <- c(min(plot_df[[column_name]]), max(plot_df[[column_name]]))
    colors <- c("blue", "red")
    col <- circlize::colorRamp2(breaks, colors)
  } else {
    stop(glue("Unknown datatype for {column_name}"))
  }

  return(col)
}

color_mappings <- lapply(names(plot_df), map_colors)
names(color_mappings) <- names(plot_df)
message(class(color_mappings))
message(names(color_mappings))
message(color_mappings)

col_ha <- ComplexHeatmap::HeatmapAnnotation(df = plot_df, col = color_mappings)
row_ha <- ComplexHeatmap::HeatmapAnnotation(df = plot_df, which = "row", show_legend = FALSE, show_annotation_name = FALSE, col = color_mappings)

# euclidean distance
pdf(snakemake@output$heatmap, width = 4 + (nsamples / 4), height = 4 + (nsamples / 4))
dst <- t(assay(vst)) %>%
  dist(diag = TRUE, upper = TRUE) %>%
  as.matrix()
p <- ComplexHeatmap::Heatmap(
  dst,
  bottom_annotation = col_ha,
  right_annotation = row_ha,
  name = "euclidean dst",
  column_title = "Euclidean Distance",
  width = nsamples * unit(5, "mm"),
  height = nsamples * unit(5, "mm")
)
print(p)

# pearson correlation
corr <- cor(assay(vst))
p <- ComplexHeatmap::Heatmap(
  corr,
  bottom_annotation = col_ha,
  right_annotation = row_ha,
  name = "pearson",
  column_title = "Pearson Correlation",
  width = nsamples * unit(5, "mm"),
  height = nsamples * unit(5, "mm")
)
print(p)
dev.off()
