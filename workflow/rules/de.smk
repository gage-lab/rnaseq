# this rule creates DESeq and SummarizedExperiment objects for DGE, DTE, and DTU analysis
rule deseq_swish:
    input:
        quant=expand(rules.salmon_quant.output.quant_tx, sample=samples["sample_name"]),
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf"),
        samplesheet=config["samples"],
    output:
        dge=f"{outdir}/dge/dge.rds",  # differential gene expression
        dte=f"{outdir}/dte/dte.rds",  # differential transcript expression
        dtu=f"{outdir}/dtu/dtu.rds",  # differential transcript usage
    log:
        f"{outdir}/deseq_swish.log",
    params:
        model=config["de"]["model"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/deseq_swish.R"


rule deseq_tetranscripts:
    input:
        quant=expand(rules.tetranscripts_quant.output, sample=samples["sample_name"]),
        txome_gtf=expand(rules.get_ref.output, file="txome.gtf"),
        rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf"),
        samplesheet=config["samples"],
    output:
        f"{outdir}/dge_te/dge_te.rds",
    log:
        log=f"{outdir}/dge_te/deseq_tetranscripts.log",
    params:
        model=config["de"]["model"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/deseq_tetrancripts.R"


rule pca_heatmap:
    input:
        rules.deseq_swish.output,
    output:
        pca_plot=f"{outdir}/{{de}}/pca.svg",
        scree_plot=f"{outdir}/{{de}}/scree.svg",
        pairs_plot=f"{outdir}/{{de}}/pairs.svg",
        dist_heatmap=f"{outdir}/{{de}}/dist_heatmap.svg",
        corr_heatmap=f"{outdir}/{{de}}/corr_heatmap.svg",
    log:
        f"{outdir}/{{de}}/pca_heatmap.log",
    params:
        model=config["de"]["model"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/pca_heatmap.R"


rule results:
    input:
        f"{outdir}/{{de}}/{{de}}.rds",
    output:
        f"{outdir}/{{de}}/{{contrast}}_results.csv",
    log:
        f"{outdir}/{{de}}/{{contrast}}_results.log",
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
        MA_plot=f"{outdir}/{{de}}/{{contrast}}_MA.svg",
        volcano_plot=f"{outdir}/{{de}}/{{contrast}}_volcano.svg",
    log:
        f"{outdir}/{{de}}/{{contrast}}_volcano_MA.log",
    params:
        LFCcutoff=config["de"]["cutoffs"]["log2FoldChange"],
        FDRcutoff=config["de"]["cutoffs"]["FDR"],
        max_overlaps=15,
        label=lambda wc: "gene_name"
        if wc.de == "dge" or "dge_te"
        else "transcript_name",
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/volcano_MA.R"


rule make_gs_df:
    output:
        f"{outdir}/dge/gsea/{{gs}}_genesets.tsv",
    log:
        f"{outdir}/dge/gsea/{{gs}}_genesets.log",
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
        results=f"{outdir}/dge/gsea/{{gs}}/{{contrast}}_results.tsv",
    log:
        f"{outdir}/dge/gsea/{{gs}}/{{contrast}}.log",
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/gsea.R"


rule ora:
    input:
        dge=expand(rules.results.output, de="dge", allow_missing=True),
    output:
        resultsUP=f"{outdir}/dge/ora/{{contrast}}_UP.tsv",
        resultsDOWN=f"{outdir}/dge/ora/{{contrast}}_DOWN.tsv",
    log:
        f"{outdir}/dge/ora/{{contrast}}.log",
    params:
        LFCcutoff=config["de"]["cutoffs"]["log2FoldChange"],
        FDRcutoff=config["de"]["cutoffs"]["FDR"],
    conda:
        "../envs/de.yaml"
    script:
        "../scripts/ora.R"
