rule tetranscripts_count:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
    output:
        "{outdir}/map_count/tetranscripts/{sample}/TEtranscripts_out.cntTable",
    conda:
        "../envs/tetranscripts.yaml"
    shadow:
        "shallow"
    log:
        "{outdir}/map_count/tetranscripts/{sample}/TEtranscripts.err",
    params:
        strandedness=get_strandedness,
        mode="multi",
    shell:
        """
        touch {log} && exec >> {log} 2>&1
        mkdir -p $(dirname {output})
        TEcount \
            -b {input.bam} \
            --GTF {input.txome_gtf} --TE {input.rmsk_gtf} \
            --mode {params.mode} \
            --stranded {params.strandedness} \
            --sortByPos --verbose 3 \
            --outdir $(dirname {output})
        """


# convert counts into TPM, make compatible with tximport
rule tetranscripts_quant:
    input:
        counts=rules.tetranscripts_count.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
    output:
        "{outdir}/map_count/tetranscripts/{sample}/TEtranscripts_out.quant",
    log:
        "{outdir}/map_count/tetranscripts/{sample}/quant.err",
    conda:
        "../envs/tetranscripts.yaml"
    script:
        "../scripts/te_quant.py"
