# Gage Lab RNA-seq pipeline

[![Tests](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml/badge.svg)](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.22.0-brightgreen.svg)](https://snakemake.github.io)

Analysis pipeline for single-end or paired-end Illumina RNA-seq experiments

## Current features

1. Download reference genome, transcriptome, gene/transcript annotation, and TE annotation
2. Generating STAR index
3. Read alignment with STAR to genome and transcriptome, flexible to single-end or paired-end reads
4. Quantification of gene and transcript abundances with Salmnon
5. Quantification of TE subfamily abundances with TEtranscripts
6. Differential gene and TE expression analysis with DESeq2 (`dge`, and `dge_te`)
7. Differential transcript expression and usage analysis with Swish (`dte` and `dtu`)
8. PCA and heatmap of gene and TE expression
9. Volcano and MA plotting of all 4 differential analyses
10. Functional enrichment analysis using Gene set enrichment analysis (GSEA) and Overrepresentation analysis (ORA, aka GO analysis)

## TODO

- [ ] fix TElocal testing (running out of memory)
- [ ] Put volcano and MA plots into one pdf per contrast per diff analysis
- [ ] Fix QC to show each run for STAR and Salmon
- [ ] get strandedness from salmon output instead of samples.tsv
- [ ] add usage and contribution docs
- [ ] Change print statements to message statements in R scripts (they get printed to the log file)
