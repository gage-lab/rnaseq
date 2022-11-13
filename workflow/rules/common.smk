import pandas as pd
import numpy as np

samples = (
    pd.read_csv(config["samples"], dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

# choose reference genome fasta and gene annotation file
if config["ref"]["region"] == "chr21":
    genome = {
        "fa": ".test/ngs-test-data/ref/genome.chr21.fa",
        "fai": ".test/ngs-test-data/ref/genome.chr21.fa.fai",
        "index": directory("resources/star_genome_chr21"),
    }
    genes = ".test/ngs-test-data/ref/annotation.chr21.gtf"
elif config["ref"]["region"] == None:
    genome = {
        "fa": "resources/genome.fa",
        "fai": "resources/genome.fa.fai",
        "index": directory("resources/star_genome"),
    }
    genes = "resources/genes.gtf"
else:
    raise ValueError("Invalid reference genome region")

if "genes" in config["ref"]:
    genes = config["ref"]["genes"]


def get_bam(wildcards):
    s = samples.loc[wildcards.sample].dropna()
    if "bam" in s.index:
        return s["bam"]
    else:
        return (
            f"{config['outdir']}/star/{wildcards.sample}/Aligned.sortedByCoord.out.bam"
        )


def get_strandedness(wildcards):
    s = samples.loc[wildcards.sample]["strandedness"]
    if pd.isnull(s) or s == "unstranded":
        return "no"
    elif s == "forward":
        return "forward"
    elif s == "reverse":
        return "reverse"
    else:
        raise ValueError("Invalid strandedness value")
