# TODO: get strandedness from salmon output instead of samples.tsv
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
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
    output:
        "{outdir}/map_count/{sample}/tetranscripts/TEtranscripts_out.cntTable",
    conda:
        "../envs/tetranscripts.yaml"
    shadow:
        "shallow"
    log:
        "{outdir}/map_count/{sample}/tetranscripts/TEtranscripts.err",
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
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
    output:
        "{outdir}/map_count/{sample}/tetranscripts/TEtranscripts_out.quant",
    log:
        "{outdir}/map_count/{sample}/tetranscripts/quant.err",
    conda:
        "../envs/tetranscripts.yaml"
    script:
        "../scripts/tetranscripts_quant.py"
