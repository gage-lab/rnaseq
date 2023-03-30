rule telocal_count:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_ind=expand(rules.get_ref.output, file="rmsk.locInd", allow_missing=True),
    output:
        "{outdir}/map_count/{sample}/telocal/TElocal_out.cntTable",
    conda:
        "../envs/telocal.yaml"
    shadow:
        "shallow"
    log:
        "{outdir}/map_count/{sample}/telocal/telocal.err",
    params:
        strandedness=get_strandedness,
        mode="multi",
    shell:
        """
        touch {log} && exec >> {log} 2>&1
        mkdir -p $(dirname {output})
        TElocal \
            -b {input.bam} \
            --GTF {input.txome_gtf} --TE {input.rmsk_ind} \
            --mode {params.mode} \
            --project $(dirname {output})/$(basename {output} -s .cntTable) \
            --stranded {params.strandedness} \
            --sortByPos --verbose 3
        """


# convert counts into TPM, make compatible with tximport
rule telocal_quant:
    input:
        counts=rules.telocal_count.output,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
    output:
        "{outdir}/map_count/{sample}/telocal/TElocal_out.quant",
    log:
        "{outdir}/map_count/{sample}/telocal/quant.err",
    conda:
        "../envs/telocal.yaml"
    script:
        "../scripts/te_quant.py"
