def get_star_input(wildcards):
    if config["trimming"]["activate"]:
        if not is_paired_end(wildcards.sample):
            return {
                "fq1": f"{wildcards.outdir}/trimmed/{wildcards.sample}_trimmed.fq.gz"
            }
        else:
            return {
                "fq1": f"{wildcards.outdir}/trimmed/{wildcards.sample}_val_1.fq.gz",
                "fq2": f"{wildcards.outdir}/trimmed/{wildcards.sample}_val_2.fq.gz",
            }
    else:
        if not is_paired_end(wildcards.sample):
            s = samples.loc[(wildcards.sample), ["fq1"]].dropna()
            return {"fq1": f"{s.fq1}"}
        else:
            s = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
            return {"fq1": f"{s.fq1}", "fq2": f"{s.fq2}"}


rule star_index:
    input:
        fasta=expand(rules.get_ref.output, file="genome.fa", allow_missing=True),
    output:
        star_index=directory("{outdir}/resources/star_index"),
    threads: 8
    log:
        "{outdir}/resources/star_index.log",
    wrapper:
        "v1.20.0/bio/star/index"


rule star_align:
    input:
        unpack(get_star_input),
        idx=rules.star_index.output,
        gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
    output:
        genome_bam="{outdir}/map_count/{sample}/star/Aligned.out.bam",
        txome_bam="{outdir}/map_count/{sample}/star/Aligned.toTranscriptome.out.bam",
        log="{outdir}/map_count/{sample}/star/Log.out",
        log_final="{outdir}/map_count/{sample}/star/Log.final.out",
    threads: 8
    log:
        "{outdir}/map_count/{sample}/star/Log.err",
    params:
        # these parameters are optimized to retain multimapping reads
        # TODO: add description of each parameter
        extra=f"""--outSAMmultNmax -1 --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --outMultimapperOrder Random --runRNGseed 777 --outSAMtype BAM Unsorted --sjdbScore 1""",
    conda:
        "../envs/star.yaml"
    script:
        # use custom script to return both genome_bam and txome_bam
        "../scripts/star.py"


rule samtools_sort:
    input:
        rules.star_align.output.genome_bam,
    output:
        "{outdir}/map_count/{sample}/star/Aligned.out.sorted.bam",
    log:
        "{outdir}/map_count/{sample}/star/samtools_sort.log",
    wrapper:
        "v1.20.0/bio/samtools/sort"


rule samtools_index:
    input:
        rules.samtools_sort.output,
    output:
        "{outdir}/map_count/{sample}/star/Aligned.out.sorted.bam.bai",
    log:
        "{outdir}/map_count/{sample}/star/samtools_index.log",
    wrapper:
        "v1.20.0/bio/samtools/index"
