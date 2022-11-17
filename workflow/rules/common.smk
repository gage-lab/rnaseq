import pandas as pd
import numpy as np

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
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

# choose STAR and TEcount parameters (STAR params adapted from ENCODE parameters https://github.com/ENCODE-DCC/rna-seq-pipeline/blob/dev/src/align.py)
config["star"] = {}
config["star"][
    "extra"
] = f"""--outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 200 \
    --outMultimapperOrder Random \
    --runRNGseed 777 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbScore 1 \
    --sjdbGTFfile {genes}"""

if config["tecount"]["mode"] == "multi":
    config["star"]["extra"] += " --outSAMmultNmax -1"
elif config["tecount"]["mode"] == "uniq":
    config["star"]["extra"] += " --outSAMmultNmax 1"
else:
    raise ValueError("Invalid value for multimappers")


def is_paired_end(sample):
    sample_units = samples.loc[sample]
    null_units = sample_units.isnull()
    paired = ~null_units["fq2"]
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), f"invalid units for sample {sample}, must be all paired end or all single end"
    return all_paired

# TODO
# 
# replicate this function with new name
def get_trim_input(wildcards):
    sample_units = samples.loc[wildcards.sample]
    sample_name = sample_units["sample_name"]
    if sample_units["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(sample_units["fq2"]):
        # single end local sample
        return [sample_units.fq1]
    else:
        # paired end local sample
        return [sample_units.fq1, sample_units.fq2]


'''
# not TODO
#
# 1. trimmed file names
# 2. paths in config
# This is likley outdated 
def get_fastqc_input(wildcards):
    if "val" in wildcards.sample:
        return (
            f"{config['outdir']}/trimmed/{{wildcards.sample}}_val_1.fq.gz",
            f"{config['outdir']}/trimmed/{{wildcards.sample}}_val_2.fq.gz",
        )
    if "trimmed" in wildcards.sample:
        return (f"{config['outdir']}/trimmed/{{wildcards.sample}}_trimmed.fq.gz",)

    sample_units = samples.loc[wildcards.sample]
    sample_name = sample_units["sample_name"]
    if sample_units["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""
    if pd.isna(sample_units["fq2"]):
        # single end local sample
        return [sample_units.fq1]
    else:
        # paired end local sample
        return [sample_units.fq1, sample_units.fq2]'''

def get_fastqc_trimmed_input(wildcards): 
    if not is_paired_end(wildcards.sample):
        return {
            "fq1": f"{config['outdir']}/trimmed/{wildcards.sample}_trimmed.fq.gz"
        }
    else:
        return {
            "fq1": f"{config['outdir']}/trimmed/{wildcards.sample}_val_1.fq.gz",
            "fq2": f"{config['outdir']}/trimmed/{wildcards.sample}_val_2.fq.gz",
        }

def get_star_input(wildcards):
    s = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if config["trimming"]["activate"]:
        if not is_paired_end(wildcards.sample):
            return {
                "fq1": f"{config['outdir']}/trimmed/{wildcards.sample}_trimmed.fq.gz"
            }
        else:
            return {
                "fq1": f"{config['outdir']}/trimmed/{wildcards.sample}_val_1.fq.gz",
                "fq2": f"{config['outdir']}/trimmed/{wildcards.sample}_val_2.fq.gz",
            }
    else:
        if not is_paired_end(wildcards.sample):
            return {"fq1": f"{s.fq1}"}
        else:
            return {"fq1": f"{s.fq1}", "fq2": f"{s.fq2}"}

def get_multiqc_input(wildcards):
    if config["trimming"]["activate"]:
        return expand(rules.fastqc_trim.output, sample=samples['sample_name'])      
    else: 
        return None 
        


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
