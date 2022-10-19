rule trim_galore_pe:
    input:
        get_trim_input,
    output:
        "results/trimmed/{sample}_val_1.fq.gz",
        "results/trimmed/{sample}_val_2.fq.gz",
    log:
        "results/trimmed/{sample}.log",
    params:
        extra="-q 20",
    conda:
        "../envs/trim_galore.yml"
    threads: 4
    shell:
        """
        trim_galore --paired \
            --gzip\
            {params.extra} \
            --cores {threads} \
            --output_dir $(dirname {output[0]}) \
            --basename {wildcards.sample} \
            {input} 2> {log}
        """


rule trim_galore_se:
    input:
        get_trim_input,
    output:
        "results/trimmed/{sample}_trimmed.fq.gz",
    log:
        "results/trimmed/{sample}.log",
    params:
        extra="-q 20 --gzip",
    conda:
        "../envs/trim_galore.yml"
    threads: 4
    shell:
        """
        trim_galore \
            --gzip \
            {params.extra} \
            --cores {threads} \
            --basename {wildcards.sample} \
            --output_dir $(dirname {output[0]}) \
            {input} 2> {log}
        """
