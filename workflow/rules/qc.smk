def get_fastqc_input(wildcards):
    # SINGLE END
    if not is_paired_end(wildcards.sample):
        if "TSOfilter" in wildcards.suffix:
            # single end tso-filtered sample
            return rules.filterTSO_se.output.fastq
        elif "trimmed" in wildcards.suffix:
            # single end trimmed sample
            return rules.trim_galore_se.output
        else:
            # single end local sample
            return samples.loc[wildcards.sample]["fq1"]
    # PAIRED ENDS
    else:
        # paired end trimmed sample
        if "TSOfilter" in wildcards.suffix:
            if "1" in wildcards.suffix:
                return rules.filterTSO_pe.output.fastq1
            elif "2" in wildcards.suffix:
                return rules.filterTSO_pe.output.fastq2
        elif "trimmed" in wildcards.suffix:
            if "1" in wildcards.suffix:
                return rules.trim_galore_pe.output[0]
            elif "2" in wildcards.suffix:
                return rules.trim_galore_pe.output[1]
        # paired end local sample
        else:
            if wildcards.suffix == "1":
                return samples.loc[wildcards.sample]["fq1"]
            elif wildcards.suffix == "2":
                return samples.loc[wildcards.sample]["fq2"]


rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="{outdir}/qc/fastqc/{sample}_{suffix}_fastqc.html",
        zip="{outdir}/qc/fastqc/{sample}_{suffix}_fastqc.zip",
    params:
        "--quiet",
    log:
        "{outdir}/qc/fastqc/{sample}_{suffix}.log",
    threads: 1
    wrapper:
        "v1.19.0/bio/fastqc"


rule rseqc_gtf2bed:
    input:
        "{outdir}/resources/txome.gtf",
    output:
        bed="{outdir}/resources/txome.bed",
        db=temp("{outdir}/resources/txome.db"),
    log:
        "{outdir}/resources/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_readdist:
    input:
        bam=rules.samtools_sort.output,
        bed=rules.rseqc_gtf2bed.output.bed,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.read_distribution.txt",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.read_distribution.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam=rules.samtools_sort.output,
        bed=rules.rseqc_gtf2bed.output.bed,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.inner_distance_freq.inner_distance.txt",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.inner_distance.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        rules.samtools_sort.output,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.GC_plot.pdf",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.readgc.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdup:
    input:
        rules.samtools_sort.output,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.DupRate_plot.pdf",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.readdup.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_junction_annotation:
    input:
        bam=rules.samtools_sort.output,
        bed=rules.rseqc_gtf2bed.output.bed,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.junctionanno.junction.bed",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.junction_annotation.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam=rules.samtools_sort.output,
        bed=rules.rseqc_gtf2bed.output.bed,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.junctionSaturation_plot.pdf",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.junction_saturation.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_infer:
    input:
        bam=rules.samtools_sort.output,
        bed=rules.rseqc_gtf2bed.output.bed,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.infer_experiment.txt",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.infer_experiment.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_stat:
    input:
        rules.samtools_sort.output,
    output:
        "{outdir}/qc/{tso_filter}/{sample}.stats.txt",
    log:
        "{outdir}/qc/{tso_filter}/{sample}.stats.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


def get_fastqc_for_multiqc(wildcards):
    result = []
    for s in samples["sample_name"]:
        if is_paired_end(s):
            suffices = ["1", "2"]
            if config["trimming"]["activate"]:
                suffices += ["trimmed1", "trimmed1"]
            if wildcards.tso_filter == "tso_filter":
                suffices += ["TSOfiltered1", "TSOfiltered2"]
        else:
            suffices = ["1"]
            if config["trimming"]["activate"]:
                suffices += ["trimmed"]
            if wildcards.tso_filter == "tso_filter":
                suffices += ["TSOfiltered"]

        result += expand(
            rules.fastqc.output,
            sample=s,
            suffix=suffices,
            allow_missing=True,
        )
    return result


def get_quant_for_multiqc(wildcards):
    result = []
    if wildcards.tso_filter == "tso_filter" or not config["filterTSOforTE"]["activate"]:
        result += expand(
            rules.telocal_count.output,
            sample=samples["sample_name"],
            allow_missing=True,
        )
        result += expand(
            ruless.tetranscripts_count.output,
            sample=samples["sample_name"],
            allow_missing=True,
        )

    if wildcards.tso_filter == "no_filter":
        result += expand(
            rules.salmon_quant.output,
            sample=samples["sample_name"],
            allow_missing=True,
        )
    return result


rule multiqc:
    input:
        get_fastqc_for_multiqc,
        get_quant_for_multiqc,
        expand(
            rules.rseqc_readdist.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        expand(
            rules.rseqc_innerdis.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        expand(
            rules.rseqc_readgc.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        expand(
            rules.rseqc_readdup.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        expand(
            rules.rseqc_junction_annotation.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        expand(
            rules.rseqc_junction_saturation.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        expand(
            rules.rseqc_infer.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        expand(
            rules.rseqc_stat.output,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
    output:
        "{outdir}/multiqc_{tso_filter}.html",
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        "{outdir}/multiqc_{tso_filter}.log",
    wrapper:
        "v1.19.0/bio/multiqc"
