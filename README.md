# Gage Lab RNA-seq pipeline

[![Tests](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml/badge.svg)](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.22.0-brightgreen.svg)](https://snakemake.github.io)

Analysis pipeline for single-end or paired-end Illumina RNA-seq experiments

## Current features

1. Download reference genome, transcriptome, gene/transcript annotation, and TE annotation
2. Generating STAR index
3. Read alignment with STAR to genome and transcriptome, flexible to single-end or paired-end reads
4. Quantification of gene and transcript abundances with Salmnon
5. Quantification of TE subfamily abundances with TEtranscripts and TE locus abundances with TEcount
6. Differential gene and TE expression analysis with DESeq2 (`dge`, `dge_te_subfamily`, `dge_te_locus`)
7. Differential transcript expression and usage analysis with Swish (`dte` and `dtu`)
8. PCA and heatmap of gene and TE expression
9. Volcano and MA plotting for each differential analysis
10. Functional enrichment analysis on differential gene expression analysis using Gene set enrichment analysis (GSEA) and Overrepresentation analysis (ORA, aka GO analysis)

## TODO

- [ ] add config parameters to specify unique or multi for TEtranscripts and TElocal
- [ ] add usage and contribution docs
- [ ] add rnaseqc
- [ ] add hyperlinks to each tool used in the pipeline
