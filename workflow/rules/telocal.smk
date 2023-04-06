rule telocal_count:
    input:
        unpack(get_te_bams),
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk=expand(rules.get_ref.output, file="rmsk.locInd", allow_missing=True),
        meta_info=rules.salmon_quant.output.meta_info,
    output:
        "{outdir}/map_count/telocal/{sample}/TElocal_out.cntTable",
    conda:
        "../envs/telocal.yaml"
    shadow:
        "shallow"
    log:
        "{outdir}/map_count/telocal/{sample}/TElocal.log",
    params:
        mode="multi",  # can be multi or uniq
    script:
        "../scripts/te_count.py"


# convert counts into TPM, make compatible with tximport
rule telocal_quant:
    input:
        counts=rules.telocal_count.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
    output:
        "{outdir}/map_count/telocal/{sample}/TElocal_out.quant",
    log:
        "{outdir}/map_count/telocal/{sample}/quant.err",
    conda:
        "../envs/telocal.yaml"
    script:
        "../scripts/te_quant.py"
