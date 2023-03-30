# this rule creates DESeq and SummarizedExperiment objects for DGE, DTE, and DTU analysis
rule deseq_swish:
    input:
        quant=expand(
            rules.salmon_quant.output.quant_tx,
            sample=samples["sample_name"],
            allow_missing=True,
        ),
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        samplesheet=config["samples"],
    output:
        dge="{outdir}/de/dge.rds",  # differential gene expression
        dte="{outdir}/de/dte.rds",  # differential transcript expression
        dtu="{outdir}/de/dtu.rds",  # differential transcript usage
    log:
        "{outdir}/de/logs/deseq_swish.log",
    params:
        model=config["de"]["model"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/deseq_swish.R"


def get_deseq_te_input(wildcards):
    if wildcards.quant_level == "subfamily":
        return expand(
            rules.tetranscripts_quant.output,
            sample=samples["sample_name"],
            allow_missing=True,
        )
    elif wildcards.quant_level == "locus":
        return expand(
            rules.telocal_quant.output,
            sample=samples["sample_name"],
            allow_missing=True,
        )


rule deseq_te:
    input:
        quant=get_deseq_te_input,
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
        samplesheet=config["samples"],
    output:
        "{outdir}/de/dge_te_{quant_level}.rds",
    log:
        log="{outdir}/de/logs/deseq_te_{quant_level}.log",
    params:
        model=config["de"]["model"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/deseq_te.R"


rule pca_heatmap:
    input:
        rules.deseq_swish.output,
    output:
        pca="{outdir}/de/{de}_pca.pdf",
        heatmap="{outdir}/de/{de}_heatmaps.pdf",
    log:
        "{outdir}/de/logs/{de}_pca_heatmap.log",
    params:
        model=config["de"]["model"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/pca_heatmap.R"


rule results:
    input:
        "{outdir}/de/{de}.rds",
    output:
        "{outdir}/de/{contrast}/{de}_results.csv",
    log:
        "{outdir}/de/{contrast}/logs/{de}_results.log",
    params:
        model=config["de"]["model"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/results.R"


rule volcano_MA:
    input:
        rules.results.output,
    output:
        MA_plot="{outdir}/de/{contrast}/{de}_MA.svg",
        volcano_plot="{outdir}/de/{contrast}/logs/{de}_volcano.svg",
    log:
        "{outdir}/de/{contrast}/logs/{de}_volcano_MA.log",
    params:
        LFCcutoff=config["de"]["cutoffs"]["log2FoldChange"],
        FDRcutoff=config["de"]["cutoffs"]["FDR"],
        max_overlaps=15,
        label=lambda wc: "gene_name"
        if wc.de == "dge" or "dge_te_subfamily" or "dge_te_locus"
        else "transcript_name",
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/volcano_MA.R"


rule make_gs_df:
    output:
        "{outdir}/de/{gs}_genesets.tsv",
    log:
        "{outdir}/de/logs/{gs}_genesets.log",
    conda:
        "../envs/de.yaml"
    wildcard_constraints:
        gs="\w+",
    script:
        "../scripts/make_gs_df.R"


rule gsea:
    input:
        dge=expand(rules.results.output, de="dge", allow_missing=True),
        gs_df=rules.make_gs_df.output,
    output:
        results="{outdir}/de/{contrast}/{gs}_gsea.tsv",
    log:
        "{outdir}/de/{contrast}/logs/{gs}_gsea.log",
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/gsea.R"


rule ora:
    input:
        dge=expand(rules.results.output, de="dge", allow_missing=True),
    output:
        resultsUP="{outdir}/de/{contrast}/up_ora.tsv",
        resultsDOWN="{outdir}/de/{contrast}/down_ora.tsv",
    log:
        "{outdir}/de/{contrast}/logs/ora.log",
    params:
        LFCcutoff=config["de"]["cutoffs"]["log2FoldChange"],
        FDRcutoff=config["de"]["cutoffs"]["FDR"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/ora.R"
