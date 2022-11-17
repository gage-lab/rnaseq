def get_fastqc_trimmed_input(wildcards):
    if not is_paired_end(wildcards.sample):
        return {"fq1": f"{config['outdir']}/trimmed/{wildcards.sample}_trimmed.fq.gz"}
    else:
        return {
            "fq1": f"{config['outdir']}/trimmed/{wildcards.sample}_val_1.fq.gz",
            "fq2": f"{config['outdir']}/trimmed/{wildcards.sample}_val_2.fq.gz",
        }


rule fastqc_trim:
    input:
        f"{config['outdir']}/trimmed/{{sample}}_{{suffix}}.fq.gz",
    output:
        html=f"{config['outdir']}/fastqc/{{sample}}_{{suffix}}_fastqc_trim.html",
        zip=f"{config['outdir']}/fastqc/{{sample}}_{{suffix}}_fastqc_trim.zip",
    params:
        "--quiet",
    log:
        f"{config['outdir']}/fastqc/{{sample}}_{{suffix}}_trim.log",
    threads: 1
    wrapper:
        "v1.19.0/bio/fastqc"


def get_fastqc_notrim_input(wildcards):
    sample_units = samples.loc[wildcards.sample]
    if not is_paired_end(wildcards.sample):
        # single end local sample
        return sample_units.fq1
    else:
        # paired end local sample
        if wildcards.suffix == "1":
            return sample_units.fq1
        elif wildcards.suffix == "2":
            return sample_units.fq2


rule fastqc_notrim:
    input:
        get_fastqc_notrim_input,
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
                rules.fastqc_notrim.output,
                sample=data["sample_name"],
                suffix=["1", "2"],
            )
            if config["trimming"]["activate"]:
                result += expand(
                    rules.fastqc_trim.output,
                    sample=data["sample_name"],
                    suffix=["val_1", "val_2"],
                )
        else:
            result += expand(
                rules.fastqc_notrim.output,
                sample=data["sample_name"],
                suffix=["1"],
            )
            if config["trimming"]["activate"]:
                result += expand(
                    rules.fastqc_trim.output,
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
