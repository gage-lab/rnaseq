rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("results/pipe/cutadapt/{sample}.{fq}.{ext}"),
    log:
        "results/pipe/cutadapt/{sample}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="results/cutadapt/{sample}_R1.fastq.gz",
        fastq2="results/cutadapt/{sample}_R2.fastq.gz",
        qc="results/cutadapt/{sample}.paired.qc.txt",
    log:
        "results/cutadapt/{sample}.log",
    params:
        others=config["trimming"]["cutadapt-pe"],
        adapters=lambda w: str(samples.loc[w.sample].loc["adapters"]),
    threads: 8
    wrapper:
        "v1.17.1/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="results/cutadapt/{sample}.fastq.gz",
        qc="results/cutadapt/{sample}.single.qc.txt",
    log:
        "results/cutadapt/{sample}.log",
    params:
        others=config["trimming"]["cutadapt-se"],
        adapters_r1=lambda w: str(samples.loc[w.sample].loc["adapters"]),
    threads: 8
    wrapper:
        "v1.17.1/bio/cutadapt/se"
