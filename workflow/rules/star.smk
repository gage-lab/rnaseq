def get_star_input(wildcards):
    if wildcards.tso_filter == "tso_filter":
        if not is_paired_end(wildcards.sample):
            return {"fq1": rules.filterTSO_se.output.fastq}
        else:
            return {
                "fq1": rules.filterTSO_pe.output.fastq1,
                "fq2": rules.filterTSO_pe.output.fastq2,
            }
    elif config["trimming"]["activate"]:
        if not is_paired_end(wildcards.sample):
            return {"fq1": rules.fastp_se.output.trimmed}
        else:
            return {
                "fq1": rules.fastp_pe.output.trimmed[0],
                "fq2": rules.fastp_pe.output.trimmed[1],
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
        gtf=rules.filter_gtf.output.gtf,
    output:
        genome_bam="{outdir}/map_count/star/{sample}/{tso_filter}/Aligned.out.bam",
        txome_bam="{outdir}/map_count/star/{sample}/{tso_filter}/Aligned.toTranscriptome.out.bam",
        log="{outdir}/map_count/star/{sample}/{tso_filter}/Log.out",
        log_final="{outdir}/map_count/star/{sample}/{tso_filter}/Log.final.out",
    threads: 8
    log:
        "{outdir}/map_count/star/{sample}/{tso_filter}/Log.err",
    params:
        # allowing for a maximum of 100 multi mapping loci and 200 anchors (used by Hammell Lab)
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
        rules.star_align.output.genome_bam.replace("bam", "sorted.bam"),
    log:
        rules.star_align.log[0].replace("Log.err", "sort.log"),
    wrapper:
        "v1.20.0/bio/samtools/sort"


rule samtools_index:
    input:
        rules.samtools_sort.output,
    output:
        rules.samtools_sort.output[0] + ".bai",
    log:
        rules.samtools_sort.log[0].replace("sort", "index"),
    wrapper:
        "v1.20.0/bio/samtools/index"


# call this function for inputs to TElocal and TEtranscripts
def get_te_bams(wildcards):
    f = "tso_filter" if config["filterTSOforTE"]["activate"] else "no_filter"
    return {
        "bam": expand(rules.samtools_sort.output, tso_filter=f, allow_missing=True),
        "bai": expand(rules.samtools_index.output, tso_filter=f, allow_missing=True),
    }
