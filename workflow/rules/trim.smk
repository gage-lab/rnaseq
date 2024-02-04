def get_trim_input(wildcards):
    sample_units = samples.loc[wildcards.sample]
    if not is_paired_end(wildcards.sample):
        # single end local sample
        return {"sample": [sample_units.fq1]}
    else:
        # paired end local sample
        return {"sample": [sample_units.fq1, sample_units.fq2]}


rule fastp_se:
    input:
        unpack(get_trim_input),
    output:
        trimmed="{outdir}/trimmed/se/{sample}.fastq",
        failed="{outdir}/trimmed/se/{sample}.failed.fastq",
        html="{outdir}/trimmed/se/{sample}.html",
        json="{outdir}/trimmed/se/{sample}.json",
    log:
        "{outdir}/trimmed/se/{sample}.log",
    params:
        adapters="",
        extra="",
    threads: 1
    wrapper:
        "v3.3.5/bio/fastp"


rule fastp_pe:
    input:
        unpack(get_trim_input),
    output:
        trimmed=[
            "{outdir}/trimmed/pe/{sample}.1.fastq",
            "{outdir}/trimmed/pe/{sample}.2.fastq",
        ],
        # Unpaired reads separately
        unpaired1="{outdir}/trimmed/pe/{sample}.u1.fastq",
        unpaired2="{outdir}/trimmed/pe/{sample}.u2.fastq",
        failed="{outdir}/trimmed/pe/{sample}.failed.fastq",
        html="{outdir}/trimmed/pe/{sample}.html",
        json="{outdir}/trimmed/pe/{sample}.json",
    log:
        "{outdir}/trimmed/pe/{sample}.log",
    params:
        adapters="",
        extra="",
    threads: 2
    wrapper:
        "v3.3.5/bio/fastp"
