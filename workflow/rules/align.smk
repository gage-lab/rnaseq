rule star:
    input:
        unpack(get_star_input),
        idx=rules.star_index.output,
    output:
        aln=f"{config['outdir']}/star/{{sample}}/Aligned.sortedByCoord.out.bam",
        log=f"{config['outdir']}/star/{{sample}}/Log.out",
        log_final=f"{config['outdir']}/star/{{sample}}/Log.final.out",
    params:
        extra=config["star"]["extra"],
    threads: 8
    wrapper:
        "v1.16.0/bio/star/align"


rule samtools_index:
    input:
        get_bam,
    output:
        f"{config['outdir']}/samtools_index/{{sample}}.bam.bai",
    log:
        f"{config['outdir']}/samtools_index/{{sample}}/samtools_index.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.16.0/bio/samtools/index"
