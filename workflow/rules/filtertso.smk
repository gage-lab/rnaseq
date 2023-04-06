# https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
# https://cutadapt.readthedocs.io/en/stable/guide.html#


def get_filterTSO_pe_input(wildcards):
    if config["trimming"]["activate"]:
        return {
            "fq1": rules.trim_galore_pe.output[0],
            "fq2": rules.trim_galore_pe.output[1],
        }
    else:
        s = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
        return {"fq1": f"{s.fq1}", "fq2": f"{s.fq2}"}


rule filterTSO_pe:
    input:
        unpack(get_filterTSO_pe_input),
    output:
        fastq1="{outdir}/tso_filtered/{sample}_1.fq.gz",
        fastq2="{outdir}/tso_filtered/{sample}_2.fq.gz",
        qc="{outdir}/tso_filtered/{sample}_qc.txt",
    params:
        adapters="-g {} -G {}".format(
            config["filterTSOforTE"]["TSO"], config["filterTSOforTE"]["TSO"]
        ),
        extra="--minimum-length 20 --overlap 10 --pair-filter=both --discard-untrimmed",
    threads: 4
    wrapper:
        "v1.25.0/bio/cutadapt/pe"


def get_filterTSO_se_input(wildcards):
    if config["trimming"]["activate"]:
        return {"fq1": rules.trim_galore_se.output}
    else:
        s = samples.loc[(wildcards.sample), ["fq1"]].dropna()
        return {"fq1": f"{s.fq1}"}


rule filterTSO_se:
    input:
        unpack(get_filterTSO_se_input),
    output:
        fastq="{outdir}/tso_filtered/{sample}_filtered.fq.gz",
        qc="{outdir}/tso_filtered/{sample}_qc.txt",
    params:
        adapters="-g {}".format(config["filterTSOforTE"]["TSO"]),
        extra="--minimum-length 20 --overlap 10 --pair-filter=both --discard-untrimmed",
    threads: 4
    wrapper:
        "v1.25.0/bio/cutadapt/se"
