# TODO:
#
# 1. fastqc
# 2. MultiQC
# 3. RNASEQ


# TODO: fastqc
#
# 1. get paths for trimmed: trim.outputs
# 2. path for untrimmed should theoritically be the same
#


rule fastqc_trim:
    input:
        get_fastqc_trimmed_in
    output:
        html=f"{config['outdir']}/fastqc/{{sample}}_fastqc.html",
        zip=f"{config['outdir']}/fastqc/{{sample}}_fastqc.zip",
    params:
        "--quiet",
    log:
        f"{config['outdir']}/fastqc/{{sample}}.log",
    threads: 1
    wrapper:
        "v1.19.0/bio/fastqc"


rule fastqc_notrim: 
    input:
        get_trim_input
    output:
        html=f"{config['outdir']}/fastqc/{{sample}}_fastqc.html",
        zip=f"{config['outdir']}/fastqc/{{sample}}_fastqc.zip",
    params:
        "--quiet",
    log:
        f"{config['outdir']}/fastqc/{{sample}}.log",
    threads: 1
    wrapper:
        "v1.19.0/bio/fastqc"

rule fastqc:
    input:
        get_fastqc_input,
    output:
        html=f"{config['outdir']}/fastqc/{{sample}}_fastqc.html",
        zip=f"{config['outdir']}/fastqc/{{sample}}_fastqc.zip",
    params:
        "--quiet",
    log:
        f"{config['outdir']}/fastqc/{{sample}}.log",
    threads: 1
    wrapper:
        "v1.19.0/bio/fastqc"

rule multiqc:
    input:
        expand(rules.fastqc_notrim.output, sample=samples['sample_name'])
        #rules.fastqc_trim.output,
        #rules.star.output,

    output:
        f"{config['outdir']}/qc/multiqc.html",
    params:
        extra="",  # Optional: extra parameters for multiqc.
    log:
        f"{config['outdir']}/qc/multiqc.log",
    wrapper:
        "v1.19.0/bio/multiqc"
