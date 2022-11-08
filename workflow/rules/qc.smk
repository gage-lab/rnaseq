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

rule fastqc: 
    input:
        #rules.trim_galore_pe.output, rules.trim_galore_se.output, get_trim_input No not this 
    output:
        f"{config['outdir']}/fastqc/{{sample}}_fastqc.html", 
        f"{config['outdir']}/fastqc/{{sample}}_fastqc.zip"
    log:
        f"{config['outdir']}/fastqc/{{sample}}.log",
    conda:
        "../envs/all_qc.yml"
    wrapper:
    
