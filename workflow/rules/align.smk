rule star:
    input:
        unpack(get_fq),
        idx=rules.star_index.output,
    output:
        aln="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="results/star/{sample}/Log.out",
        log_final="results/star/{sample}/Log.final.out",
    log:
        "results/star/{sample}/Log.err",
    params:
        extra="--outSAMtype BAM SortedByCoordinate --sjdbGTFfile {}".format(
            rules.get_gencode.output
        ),
    threads: 8
    wrapper:
        "v1.16.0/bio/star/align"


rule samtools_index:
    input:
        rules.star.output.aln,
    output:
        f"{rules.star.output.aln}.bai",
    log:
        "results/star/{sample}/samtools_index.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.16.0/bio/samtools/index"
