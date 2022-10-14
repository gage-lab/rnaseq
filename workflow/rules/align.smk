rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "results/star/{sample}.log",
    params:
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {}".format(
            "resources/genome.gtf"
        ),
    threads: 8
    wrapper:
        "v1.16.0/bio/star/align"
