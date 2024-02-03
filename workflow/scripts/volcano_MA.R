#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Created on: Dec 29, 2022

con <- file(snakemake@log[[1]], "w")
sink(file = con, type = "message")

suppressPackageStartupMessages({
  library(ggplot2)
  library(glue)
  library(readr)
  library(ggrepel)
})

ggplot2::theme_set(ggplot2::theme_bw())

res <- readr::read_csv(snakemake@input[[1]])
FDRcutoff <- snakemake@config$de$cutoffs$FDR
LFCcutoff <- snakemake@config$de$cutoffs$log2FoldChange
# save.image(glue("{snakemake@wildcards[['de']]}_plot.RData")) # for debugging

n_down <- nrow(res[res$padj < FDRcutoff & res$log2FoldChange > LFCcutoff, ])
n_up <- nrow(res[res$padj < FDRcutoff & res$log2FoldChange < -LFCcutoff, ])

# MA plot
# label significant features
res$sig <- res$padj < FDRcutoff
res$sig[is.na(res$sig)] <- FALSE
pdf(snakemake@output[[1]], width = 10, height = 10)
p <- ggplot() +
  geom_point(data = res[res$sig == F, ], aes(x = log10(baseMean), y = log2FoldChange, color = sig), alpha = 0.5) +
  geom_point(data = res[res$sig == T, ], aes(x = log10(baseMean), y = log2FoldChange, color = sig), alpha = 0.5) +
  geom_text_repel(
    data = res[res$sig == T, ],
    aes(x = log10(baseMean), y = log2FoldChange, label = .data[[snakemake@params[["label"]]]]),
    max.overlaps = snakemake@params["max_overlaps"][[1]]
  ) +
  scale_color_manual(
    values = c("FALSE" = "gray", "TRUE" = "darkred"),
    labels = c(bquote("P"[adj] ~ ">" ~ .(FDRcutoff)), bquote("P"[adj] ~ "<" ~ .(FDRcutoff)))
  ) +
  labs(
    caption = glue("{n_down} features downregulated; {n_up} features upregulated"),
    color = NULL,
    x = bquote(~ Log[10] ~ "Mean"),
    y = bquote(~ Log[2] ~ "FoldChange")
  ) +
  theme(plot.title = element_text(face = "bold"))

print(p)


# volcano plot
res <- res[!is.na(res$padj), ] # remove NAs

# group points for coloring
# TODO: use colorblind-friendly colors
res$color <- character(length = nrow(res))
res$color <- "gray"
res$color[res$log2FoldChange > LFCcutoff] <- "darkred"
res$color[res$log2FoldChange < -LFCcutoff] <- "darkblue"
res$color[res$log2FoldChange > LFCcutoff & res$padj < FDRcutoff] <- "red"
res$color[res$log2FoldChange < -LFCcutoff & res$padj < FDRcutoff] <- "blue"

p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(fill = color, size = log10(baseMean)), color = "black", alpha = 0.5, pch = 21) +
  scale_fill_manual(
    values = c(
      "darkred" = "darkred",
      "darkblue" = "darkblue",
      "red" = "red",
      "blue" = "blue",
      "gray" = "gray"
    ),
    breaks = c(
      "red",
      "blue",
      "darkred",
      "darkblue"
    ),
    labels = c(
      bquote(~ Log[2] ~ "FC" ~ ">" ~ .(LFCcutoff) ~ "&" ~ "P"[adj] ~ "<" ~ .(FDRcutoff)),
      bquote(~ Log[2] ~ "FC" ~ "<" ~ -.(LFCcutoff) ~ "&" ~ "P"[adj] ~ "<" ~ .(FDRcutoff)),
      bquote(~ Log[2] ~ "FC" ~ ">" ~ .(LFCcutoff)),
      bquote(~ Log[2] ~ "FC" ~ "<" ~ -.(LFCcutoff))
    )
  ) +
  geom_text_repel(
    data = res[res$padj < FDRcutoff & abs(res$log2FoldChange) > LFCcutoff, ],
    aes(x = log2FoldChange, y = -log10(padj), label = .data[[snakemake@params[["label"]]]]),
    max.overlaps = snakemake@params["max_overlaps"][[1]],
    point_size = 0.1,
  ) +
  labs(
    caption = glue("{n_down} features downregulated; {n_up} features upregulated"),
    size = bquote(~ Log[10] ~ "Mean"),
    color = NULL,
    x = bquote(~ Log[2] ~ "FoldChange"),
    y = bquote(~ -Log[10] ~ "P"[adj]),
  ) +
  theme(plot.title = element_text(face = "bold"))

print(p)
dev.off()
