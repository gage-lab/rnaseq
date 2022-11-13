rule get_genome:
    output:
        genome["fa"],
    log:
        "resources/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "v1.17.1/bio/reference/ensembl-sequence"


# TODO: get GENCODE insteaed of Ensembl
rule get_genes:
    output:
        genes,
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "resources/get_genes.log",
    wrapper:
        "v1.17.1/bio/reference/ensembl-annotation"


rule star_index:
    input:
        fasta=rules.get_genome.output,
    output:
        genome["index"],
    threads: 8
    params:
        extra="",  # optional parameters
    log:
        "resources/star_index.log",
    wrapper:
        "v1.16.0/bio/star/index"


rule get_rmsk:
    output:
        "resources/rmsk.gtf",
    params:
        species=config["ref"]["species"],
        source=config["ref"]["source"],
        build=config["ref"]["build"],
    conda:
        "../envs/download.yml"
    log:
        "resources/get_rmsk.log",
    shell:
        """
        # TODO: update this to get most recent repeatmasker runs
        if [ {params.build} == "GRCh38" ] || [ {params.build} == "GRCm38" ]; then
            URL="https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/{params.build}_{params.source}_rmsk_TE.gtf.gz"
        else 
            echo "Species and build combination not implemented"; exit 1
        fi
        wget -O- $URL | gzip -dc > {output} 2> {log}
        """
