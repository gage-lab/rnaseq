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
        trimmed="{outdir}/trimmed/{sample}.fastq",
        failed="{outdir}/trimmed/{sample}.failed.fastq",
        html="{outdir}/trimmed/{sample}.html",
        json="{outdir}/trimmed/{sample}.json",
    log:
        "{outdir}/trimmed/{sample}.log",
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
            "{outdir}/trimmed/{sample}.1.fastq",
            "{outdir}/trimmed/{sample}.2.fastq",
        ],
        # Unpaired reads separately
        unpaired1="{outdir}/trimmed/{sample}.u1.fastq",
        unpaired2="{outdir}/trimmed/{sample}.u2.fastq",
        failed="{outdir}/trimmed/{sample}.failed.fastq",
        html="{outdir}/trimmed/{sample}.html",
        json="{outdir}/trimmed/{sample}.json",
    log:
        "{outdir}/trimmed/{sample}.log",
    params:
        adapters="",
        extra="",
    threads: 2
    wrapper:
        "v3.3.5/bio/fastp"
