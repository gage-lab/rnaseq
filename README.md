# Gage Lab RNA-seq pipeline

[![Tests](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml/badge.svg)](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.16.0-brightgreen.svg)](https://snakemake.github.io)

Alignment and quantification pipeline for single-end or paired-end Illumina RNA-seq reads

## Current features

1. Downloading reference genome, STAR index, and GENCODE annotation
2. Read alignment with STAR, flexible to single-end or paired-end reads

## Coming soon

1. handling of single-cell RNA-seq data
2. FASTQC, RNASEQC, and MultiQC reports
3. Downloading data from SRA/ENA
4. Differential expression analysis

## Development tips

Please run the following before merging to the main branch

```bash
# enforce code style
snakefmt .

# run lint checks
snakemake --lint

# test pipeline
snakemake \
   all \
   --cores 3 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --configfile .test/config/config.yaml 
```