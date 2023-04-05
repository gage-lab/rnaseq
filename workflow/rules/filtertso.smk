# https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
# https://cutadapt.readthedocs.io/en/stable/guide.html#


rule filterTSO_pe:
    input:
        rules.trim_galore_pe.output,
    output:
        fastq1="{outdir}/filtered/{sample}_1.fastq",
        fastq2="{outdir}/filtered/{sample}_2.fastq",
        qc="{outdir}/filtered/{sample}_qc.txt",
    params:
        adapters="-g {} -G {}".format(config["filterTSOforTEquant"]["sequence"]),
        extra="--minimum-length 20 --overlap 10 --pair-filter=both --discard-untrimmed",
    threads: 4
    wrapper:
        "v1.25.0/bio/cutadapt/pe"


rule filterTSO_se:
    input:
        rules.trim_galore_se.output,
    output:
        fastq="{outdir}/filtered/{sample}.fastq",
        qc="{outdir}/filtered/{sample}_qc.txt",
    params:
        adapters="-g {}".format(config["filterTSOforTEquant"]["sequence"]),
        extra="--minimum-length 20 --overlap 10 --pair-filter=both --discard-untrimmed",
    threads: 4
    wrapper:
        "v1.25.0/bio/cutadapt/se"
