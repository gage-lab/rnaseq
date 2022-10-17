rule tetranscripts_count:
    input:
        bam=rules.star.output.aln,
        bai=rules.samtools_index.output,
        gencode=rules.get_gencode.output,
        rmsk=rules.get_rmsk.output,
    output:
        directory("results/TEcount/{sample}"),
    conda:
        "../envs/tetranscripts.yml"
    shadow:
        "shallow"
    log:
        "results/TEcount/{sample}.log",
    params:
        mode=config["tecount"]["mode"],
        strandedness=get_strandedness,
    shell:
        """
        mkdir -p {output}
        TEcount \
            -b {input.bam} \
            --GTF {input.gencode} --TE {input.rmsk} \
            --mode {params.mode} \
            --stranded {params.strandedness} \
            --sortByPos --verbose 3 \
            --outdir {output} 2> {log}   
        """
