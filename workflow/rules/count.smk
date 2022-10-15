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
    shadow: "shallow"
    log:
        "results/TEcount/{sample}/Log.out",
    shell:
        """
        mkdir -p {output}
        TEcount --sortByPos --verbose 3 \
            -b {input.bam} --GTF {input.gencode} --TE {input.rmsk} \
            --outdir {output} 2> {log}   
        """
