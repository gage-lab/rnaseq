rule tetranscripts_count:
    input:
        bam=rules.star.output.aln,
        bai=rules.samtools_index.output,
        genes=rules.get_genes.output,
        rmsk=rules.get_rmsk.output,
    output:
        f"{config['outdir']}/TEcount/{{sample}}/TEtranscripts_out.cntTable",
    conda:
        "../envs/tetranscripts.yml"
    shadow:
        "shallow"
    log:
        f"{config['outdir']}/TEcount/{{sample}}/Log.err",
    params:
        mode=config["tecount"]["mode"],
        strandedness=get_strandedness,
    shell:
        """
        mkdir -p $(dirname {output})
        TEcount \
            -b {input.bam} \
            --GTF {input.genes} --TE {input.rmsk} \
            --mode {params.mode} \
            --stranded {params.strandedness} \
            --sortByPos --verbose 3 \
            --outdir $(dirname {output}) 2> {log}   
        """


rule aggregate_counts:
    input:
        expand(rules.tetranscripts_count.output, sample=samples["sample_name"]),
    output:
        f"{config['outdir']}/TEtranscripts_out.cntTable",
    log:
        f"{config['outdir']}/TEtranscripts_out.cntTable.log",
    run:
        for i, f in enumerate(input):
            df = pd.read_csv(
                f, sep="\t", index_col=0, header=0, names=[samples["sample_name"][i]]
            )
            df_all = df if i == 0 else df_all.add(df, fill_value=0)
        with open(output[0], "w") as f:
            df_all.to_csv(f, sep="\t", index_label = "id")

rule prep_deseq:
    input:
        counts = rules.aggregate_counts.output,
        rmsk = rules.get_rmsk.output,
        genes = rules.get_genes.output,
        coldata = config["samples"],
    output:
        f"{config['outdir']}/dds.rds",
    conda: 
        "../envs/deseq.yml"
    script:
        "../scripts/prep_deseq.R"