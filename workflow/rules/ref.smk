def get_ref_input(wildcards):
    if config["ref"][wildcards.file].startswith("ftp"):
        return FTP.remote(config["ref"][wildcards.file], static=True)
    elif config["ref"][wildcards.file].startswith("http"):
        return HTTP.remote(config["ref"][wildcards.file], static=True)
    elif os.path.exists(config["ref"][wildcards.file]):
        return config["ref"][wildcards.file]
    else:
        return AUTO.remote(config["ref"][wildcards.file])


rule get_ref:
    input:
        get_ref_input,
    output:
        "{outdir}/resources/{file}",
    log:
        "{outdir}/resources/get_{file}.log",
    conda:
        "../envs/get_ref.yaml"
    shell:
        """
        touch {log} && exec 2>&1 1>>{log}

        if [[ {input} == *"gz" ]]; then
            gzip -dc {input} > {output}
        else
            rsync {input} {output}
        fi
        """


rule rmsk_genes:
    input:
        genome_fa=expand(rules.get_ref.output, file="genome.fa", allow_missing=True),
        txome_fa=expand(rules.get_ref.output, file="txome.fa", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
    output:
        fa="{outdir}/resources/rmsk_genes.fa",
        gtf="{outdir}/resources/rmsk_genes.gtf",
    log:
        "{outdir}/resources/rmsk_genes.log",
    conda:
        "../envs/get_ref.yaml"
    shell:
        """
        touch {log} && exec 2>&1 1>>{log}

        gffread -g {input.genome_fa} -w $TMP/rmsk.fa {input.rmsk_gtf}
        cat $TMP/rmsk.fa {input.txome_fa} > {output.fa}
        cat {input.rmsk_gtf} {input.txome_gtf} > {output.gtf}
        """
