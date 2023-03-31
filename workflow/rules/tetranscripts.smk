rule tetranscripts_count:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
        meta_info=rules.salmon_quant.output.meta_info,
    output:
        "{outdir}/map_count/tetranscripts/{sample}/TEtranscripts_out.cntTable",
    conda:
        "../envs/tetranscripts.yaml"
    shadow:
        "shallow"
    log:
        "{outdir}/map_count/tetranscripts/{sample}/TEtranscripts.log",
    params:
        mode="multi",  # can be multi or uniq
    script:
        "../scripts/te_count.py"


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
