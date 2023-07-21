rule get_l1em:
    output:
        directory("{outdir}/resources/L1EM"),
    shell:
        "git clone https://github.com/FenyoLab/L1EM {output}"


rule bwa_index:
    input:
        expand(rules.get_ref.output, file="genome.fa", allow_missing=True),
    output:
        idx=multiext(
            "{outdir}/resources/genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "{outdir}/resources/bwa_index.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.28.0/bio/bwa/index"


rule build_l1em_ref:
    input:
        genome=expand(rules.get_ref.output, file="genome.fa", allow_missing=True),
        index=rules.bwa_index.output,
        l1em=rules.get_l1em.output,
    output:
        "{outdir}/resources/build_l1em_ref.log",
    params:
        species=config["species"],
    shell:
        """
        cd {input.l1em}

        # if species is human
        if [[ {params.species} == "human" ]]; then
            bash generate_L1EM_fasta_and_index.sh {input.genome} 2> {output}
        elif [[ {params.species} == "mouse" ]]; then
            bash generate_mm39_L1EM_fasta_and_index.sh {input.genome} 2> {output}
        fi
        """


rule l1em:
    input:
        bam=expand(
            rules.samtools_sort.output, tso_filter="no_filter", allow_missing=True
        ),
        bai=expand(
            rules.samtools_index.output, tso_filter="no_filter", allow_missing=True
        ),
        l1em=rules.get_l1em.output,
        index=rules.bwa_index.output,
        ref=expand(rules.get_ref.output, file="genome.fa", allow_missing=True),
        l1em_ref=rules.build_l1em_ref.output,
    output:
        full_counts="{outdir}/map_count/L1EM/{sample}/full_counts.txt",
        l1hs_counts="{outdir}/map_count/L1EM/{sample}/l1hs_transcript_counts.txt",
        filter_fpm="{outdir}/map_count/L1EM/{sample}/filter_L1HS_FPM.txt",
    threads: 2
    conda:
        "../envs/L1EM.yaml"
    params:
        species=config["species"],
    log:
        "{outdir}/map_count/L1EM/{sample}/l1em.log",
    shell:
        """
        touch {log} && exec >> {log} 2>&1
        outdir=$(dirname {output.full_counts})
        mkdir -p $outdir && cd $outdir

        if [[ {params.species} == "human" ]]; then
            bash -e {input.l1em}/run_L1EM.sh {input.bam} {input.l1em} {input.ref}
        elif [[ {params.species} == "mouse" ]]; then
            bash -e {input.l1em}/run_L1EM_mm39.sh {input.bam} {input.l1em} {input.ref}
        fi

        """
