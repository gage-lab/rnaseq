rule rnaseq_pipeline:
    input:
        input=config["samples"],
        fasta=rules.get_genome.output,
        gtf=rules.get_genes.output,
        star_index=rules.star_index.output,
    output:
        multiqc="{outdir}/multiqc/star_rsem/{title}_multiqc_report.html".format(
            outdir=config["outdir"], title=config["nextflow"]["multiqc_title"]
        ),
    log:
        "{outdir}/nfcore_rnaseq.log".format(outdir=config["outdir"]),
    handover: True
    threads: 8
    params:
        **config["nextflow"],
    wrapper:
        "file://workflow/wrappers/nextflow"


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
