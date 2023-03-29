rule salmon_quant:
    input:
        bam=rules.star_align.output.txome_bam,
        txome=expand(rules.get_ref.output, file="txome.fa"),
        gtf=expand(rules.get_ref.output, file="txome.gtf"),
    output:
        quant_tx=f"{outdir}/salmon_quant/{{sample}}/quant.sf",
        quant_ge=f"{outdir}/salmon_quant/{{sample}}/quant.genes.sf",
    log:
        f"{outdir}/salmon_quant/{{sample}}/logs/salmon_quant.log",
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


def get_strandedness(wildcards):
    s = samples.loc[wildcards.sample]["strandedness"]
    if pd.isnull(s) or s == "unstranded":
        return "no"
    elif s == "forward":
        return "forward"
    elif s == "reverse":
        return "reverse"
    else:
        raise ValueError("Invalid strandedness value")


rule tetranscripts_count:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf"),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf"),
    output:
        f"{outdir}/TEcount/{{sample}}/TEtranscripts_out.cntTable",
    conda:
        "../envs/tetranscripts.yaml"
    shadow:
        "shallow"
    log:
        f"{outdir}/TEcount/{{sample}}/count.err",
    params:
        strandedness=get_strandedness,
        mode="multi",
    shell:
        """
        mkdir -p $(dirname {output})
        TEcount \
            -b {input.bam} \
            --GTF {input.txome_gtf} --TE {input.rmsk_gtf} \
            --mode {params.mode} \
            --stranded {params.strandedness} \
            --sortByPos --verbose 3 \
            --outdir $(dirname {output}) 2> {log}   
        """


rule tetranscripts_quant:
    input:
        counts=rules.tetranscripts_count.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf"),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf"),
    output:
        f"{outdir}/TEcount/{{sample}}/TEtranscripts_out.quant",
    log:
        f"{outdir}/TEcount/{{sample}}/quant.err",
    conda:
        "../envs/tetranscripts.yaml"
    script:
        "../scripts/tetranscripts_quant.py"
