import pandas as pd

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
    }
    genes = ".test/ngs-test-data/ref/annotation.chr21.gtf"
elif config["ref"]["region"] == None:
    genome = {"fa": "resources/genome.fa", "fai": "resources/genome.fa.fai"}
    genes = "resources/gencode.gtf"
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
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def get_fq(wildcards):
    s = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if not is_paired_end(wildcards.sample):
        return {"fq1": f"{s.fq1}"}
    else:
        return {"fq1": f"{s.fq1}", "fq2": f"{s.fq2}"}
