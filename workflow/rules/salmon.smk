rule salmon_quant:
    input:
        bam=rules.star_align.output.txome_bam,
        txome=expand(rules.get_ref.output, file="txome.fa", allow_missing=True),
        gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
    output:
        quant_tx="{outdir}/map_count/{sample}/quant.sf",
        quant_ge="{outdir}/map_count/{sample}/quant.genes.sf",
    log:
        "{outdir}/map_count/{sample}/logs/salmon_quant.log",
    params:
        libtype="A",
        numBootstraps=config["salmon"]["numBootstraps"],
        extra=config["salmon"]["extra"],
    threads: 2
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant \
            --threads {threads} \
            --geneMap {input.gtf} \
            --libType {params.libtype} \
            --numBootstraps {params.numBootstraps} \
            -t {input.txome} \
            -a {input.bam} \
            -o $(dirname {output.quant_tx}) \
        {params.extra} &> /dev/null
        """
