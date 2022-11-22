def get_fastqc_input(wildcards):
    # Adjust wildcards if snakemake mistakes with underscores
    if "val" in wildcards.sample:
        wildcards.suffix = wildcards.sample[-4:] + wildcards.suffix
        wildcards.sample = wildcards.sample[:-4]

    # SINGLE END
    if not is_paired_end(wildcards.sample):
        if "trimmed" in wildcards.suffix:
            # single end trimmed sample
            return f"{config['outdir']}/trimmed/{wildcards.sample}_trimmed.fq.gz"
        else:
            # single end local sample
            return samples.loc[wildcards.sample]["fq1"]
    # PAIRED ENDS
    else:
        # paired end local sample
        if "val" in wildcards.suffix:
            if "1" in wildcards.suffix:
                return f"{config['outdir']}/trimmed/{wildcards.sample}_val_1.fq.gz"

            elif "2" in wildcards.suffix:
                return f"{config['outdir']}/trimmed/{wildcards.sample}_val_2.fq.gz"
        # paired end local sample
        else:
            if wildcards.suffix == "1":
                return samples.loc[wildcards.sample]["fq1"]
            elif wildcards.suffix == "2":
                return samples.loc[wildcards.sample]["fq2"]


rule fastqc:
    input:
        get_fastqc_input,
    output:
        html=f"{config['outdir']}/fastqc/{{sample}}_{{suffix}}_fastqc.html",
        zip=f"{config['outdir']}/fastqc/{{sample}}_{{suffix}}_fastqc.zip",
    params:
        "--quiet",
    log:
        f"{config['outdir']}/fastqc/{{sample}}_{{suffix}}.log",
    threads: 1
    wrapper:
        "v1.19.0/bio/fastqc"


def get_multiqc_input(wildcards):
    result = []
    # itr over samples, is it paired or single?, adjust suffix based off of trim name
    for index, data in samples.iterrows():
        if is_paired_end(data["sample_name"]):
            result += expand(
                rules.fastqc.output,
                sample=data["sample_name"],
                suffix=["1", "2"],
            )
            if config["trimming"]["activate"]:
                result += expand(
                    rules.fastqc.output,
                    sample=data["sample_name"],
                    suffix=["val_1", "val_2"],
                )
        else:
            result += expand(
                rules.fastqc.output,
                sample=data["sample_name"],
                suffix=["1"],
            )
            if config["trimming"]["activate"]:
                result += expand(
                    rules.fastqc.output,
                    sample=data["sample_name"],
                    suffix=["trimmed"],
                )
    return result


rule multiqc:
    input:
        get_multiqc_input,
        expand(rules.star.output, sample=samples["sample_name"]),
    output:
        f"{config['outdir']}/multi/multiqc.html",
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        f"{config['outdir']}/multi/multiqc.log",
    wrapper:
        "v1.19.0/bio/multiqc"
