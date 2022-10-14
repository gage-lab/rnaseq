rule get_genome:
	output:
		fa = "resources/genome.fa",
		fai = "resources/genome.fa.fai",
		dictionary = "resources/genome.dict",
		xml = "resources/GenomeSize.xml",
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

rule get_star_index:
	output:
		directory("resources/star_genome"),
	log:
		"resources/star_genome/get_star_index.log",
	params:
		species=config["ref"]["species"],
		build=config["ref"]["build"],
	conda: 
		"../envs/download.yml"
	shell:
		"""
		URL="s3://ngi-igenomes/igenomes/{params.species}/UCSC/{params.build}/Sequence/STARIndex/"
		aws s3 --no-sign-request --region eu-west-1 sync $URL {output} 2> {log}
		"""

rule get_gencode:
	output:
		"resources/gencode.gtf",
	params:
		species=config["ref"]["species"],
		build=config["ref"]["build"],
	conda: "../envs/download.yml"
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

