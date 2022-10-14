rule tetranscripts_count:
    input:
        bam=rules.star.output.aln,
        bai=rules.samtools_index.output,
        gencode=rules.get_gencode.output,
        rmsk=rules.get_rmsk.output,
    output:
        directory("results/TEcount/{sample}"),
    container:
        "docker://mhammelllab/tetranscripts:2.2.3"
    log:
        "results/TEcount/{sample}.log",
    shell:
        """
        mkdir -p {output}
        TEcount -b {input.bam} --GTF {input.gencode} --TE {input.rmsk} --sortByPos --outdir {output} 2> {log}   
        """
