rule get_genome:
    output:
        **genome,
    log:
        "resources/get-genome.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
    conda:
        "../envs/download.yml"
    shell:
        """
        URL=s3://ngi-igenomes/igenomes/{params.species}/UCSC/{params.build}/Sequence/WholeGenomeFasta/
        aws s3 --no-sign-request --region eu-west-1 sync $URL ./resources 2> {log}
        """


rule star_index:
    input:
        fasta=rules.get_genome.output.fa,
    output:
        directory("resources/star_genome"),
    threads: 8
    params:
        extra="",
    log:
        "resources/star_index.log",
    wrapper:
        "v1.16.0/bio/star/index"


rule get_gencode:
    output:
        genes,
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
    conda:
        "../envs/download.yml"
    log:
        "resources/get_gencode.log",
    shell:
        """
        if [ {params.species} == "Homo_sapiens" ] && [ {params.build} == "hg38" ]; then
            URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
        elif [ {params.species} == "Mus_musculus" ] && [ {params.build} == "mm39" ]; then
            URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
        else 
            echo "Species and build combination not implemented"; exit 1
        fi
        wget -O- $URL | gzip -dc > {output} 2> {log}
        """


rule get_rmsk:
    output:
        "resources/rmsk.gtf",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
    conda:
        "../envs/download.yml"
    log:
        "resources/get_rmsk.log",
    shell:
        """
        # TODO: update this to get most recent repeatmasker runs
        if [ {params.species} == "Homo_sapiens" ] && [ {params.build} == "hg38" ]; then
            URL="https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_GENCODE_rmsk_TE.gtf.gz"
        elif [ {params.species} == "Mus_musculus" ] && [ {params.build} == "mm39" ]; then
            URL="https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm39_GENCODE_rmsk_TE.gtf.gz"
        else 
            echo "Species and build combination not implemented"; exit 1
        fi
        wget -O- $URL | gzip -dc > {output} 2> {log}
        """
