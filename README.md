# Gage Lab RNA-seq pipeline

![Tests](https://github.com/gage-lab/rnaseq/actions/workflows/main.yml/badge.svg)

Alignment and quantification pipeline for single-end or paired-end Illumina RNA-seq reads

## Current features

1. Downloading reference genome, STAR index, and GENCODE annotation
2. Read alignment with STAR, flexible to single-end or paired-end reads

## Coming soon

1. tetranscripts for quantification
2. handling of STAR parameters depending on strandedness of sequencing data

## Testing

```bash
snakemake \
   all \
   --cores 1 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --configfile .test/config/config.yaml 
```