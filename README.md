# Gage Lab RNA-seq pipeline

![Tests](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml/badge.svg)

Alignment and quantification pipeline for single-end or paired-end Illumina RNA-seq reads

## Current features

1. Downloading reference genome, STAR index, and GENCODE annotation
2. Read alignment with STAR, flexible to single-end or paired-end reads

## Coming soon

1. handling of STAR parameters depending on strandedness of sequencing data
2. handling of single-cell RNA-seq data

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
   --cores 1 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --configfile .test/config/config.yaml 
```