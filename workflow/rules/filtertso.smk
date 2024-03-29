# https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
# https://cutadapt.readthedocs.io/en/stable/guide.html#


def get_filterTSO_pe_input(wildcards):
    if config["trimming"]["activate"]:
        return {
            "fq1": rules.fastp_pe.output.trimmed[0],
            "fq2": rules.fastp_pe.output.trimmed[1],
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
    params:
        tso=config["filterTSOforTE"]["TSO"],
    conda:
        "../envs/cutadapt.yaml"
    log:
        "{outdir}/tso_filtered/{sample}.log",
    threads: 4
    shell:
        """
        tmp=$(mktemp -d)

        # get read pairs where at least one mate has 10bp of the TSO
        cutadapt -j {threads} -g {params.tso} -G {params.tso} \
            --pair-filter=both --discard-untrimmed \
            -o $tmp/1.fq.gz -p $tmp/2.fq.gz {input.fq1} {input.fq2} > {log}

        # remove read pairs that are too short
        cutadapt -j {threads} --minimum-length 20 --pair-filter=any \
            -o {output.fastq1} -p {output.fastq2} $tmp/1.fq.gz $tmp/2.fq.gz >> {log}
        """


def get_filterTSO_se_input(wildcards):
    if config["trimming"]["activate"]:
        return rules.fastp_se.output.trimmed
    else:
        s = samples.loc[(wildcards.sample), ["fq1"]].dropna()
        return s.fq1


rule filterTSO_se:
    input:
        get_filterTSO_se_input,
    output:
        fastq="{outdir}/tso_filtered/{sample}_filtered.fq.gz",
    params:
        tso=config["filterTSOforTE"]["TSO"],
    conda:
        "../envs/cutadapt.yaml"
    log:
        "{outdir}/tso_filtered/{sample}.log",
    threads: 4
    shell:
        """
        # get reads with TSO
        cutadapt -j {threads} -g {params.tso} --minimum-length 20 \
            --discard-untrimmed \
            -o {output.fastq} {input} > {log}
        """
