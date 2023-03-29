def get_trim_input(wildcards):
    sample_units = samples.loc[wildcards.sample]
    if not is_paired_end(wildcards.sample):
        # single end local sample
        return [sample_units.fq1]
    else:
        # paired end local sample
        return [sample_units.fq1, sample_units.fq2]


rule trim_galore_pe:
    input:
        get_trim_input,
    output:
        "{outdir}/trimmed/{sample}_val_1.fq.gz",
        "{outdir}/trimmed/{sample}_val_2.fq.gz",
    log:
        "{outdir}/trimmed/{sample}.log",
    params:
        extra="-q 20",  # -q 20 is quality cutoff
    threads: 4
    conda:
        "../envs/trim_galore.yaml"
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
        "{outdir}/trimmed/{sample}_trimmed.fq.gz",
    log:
        "{outdir}/trimmed/{sample}.log",
    params:
        extra="-q 20",  # -q 20 is quality cutoff
    threads: 4
    conda:
        "../envs/trim_galore.yaml"
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
