rule rnaseq_pipeline:
    input:
        input=config["rnaseq_samples"],
        fasta=rules.get_genome.output,
        gtf=rules.get_genes.output,
        star_index=rules.star_index.output,
    output:
        star=directory("star_rsem")
		multiqc="{outdir}/multiqc/star_rsem/{title}_multiqc_report.html".format(
			outdir=config["outdir"],
            title=config["nextflow"]["multiqc_title"]
        ),
    log:
        "results/nfcore_rnaseq.log",
    handover: True
    params:
        **config["nextflow"],
    shell:
	    "nextflow run -profile {params.profile} -params-file {params.params_file} "


rule cleanup_nextflow:
    input:
        rules.rnaseq_pipeline.output,
    output:
        "results/.nextflow.log",
    shell:
        """
        rm -rf work
        mv .nextflow.log {output}
        """