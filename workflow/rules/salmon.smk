rule salmon_quant:
    input:
        bam=expand(
            rules.star_align.output.txome_bam,
            tso_filter="no_filter",
            allow_missing=True,
        ),
        txome=rules.filter_gtf.output.fa,
        gtf=rules.filter_gtf.output.gtf,
    output:
        meta_info="{outdir}/map_count/salmon/{sample}/aux_info/meta_info.json",
        quant_tx="{outdir}/map_count/salmon/{sample}/quant.sf",
        quant_ge="{outdir}/map_count/salmon/{sample}/quant.genes.sf",
    log:
        "{outdir}/map_count/salmon/{sample}/log.orr",
    params:
        libtype="A",
        numBootstraps=config["salmon"]["numBootstraps"],
        extra=config["salmon"]["extra"],
    threads: 2
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        touch {log} && exec >> {log} 2>&1
        salmon quant \
            --threads {threads} \
            --geneMap {input.gtf} \
            --libType {params.libtype} \
            --numBootstraps {params.numBootstraps} \
            -t {input.txome} \
            -a {input.bam} \
            -o $(dirname {output.quant_tx}) \
        {params.extra}
        """
