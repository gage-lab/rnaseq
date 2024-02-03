storage ftp:
	provider="ftp"

storage http:
	provider="http"

def get_ref_input(wildcards):
	if config["ref"][wildcards.file].startswith("ftp"):
		return storage.ftp(
			config["ref"][wildcards.file],
			static=True,
			keep_local=True,
			immediate_close=True,
		)
	elif config["ref"][wildcards.file].startswith("http"):
		return storage.http(config["ref"][wildcards.file], static=True, keep_local=True)
	elif os.path.exists(config["ref"][wildcards.file]):
		return config["ref"][wildcards.file]
	else:
		raise ValueError(f"Could not find {wildcards.file}")


rule get_ref:
	input:
		get_ref_input,
	output:
		"{outdir}/resources/{file}",
	log:
		"{outdir}/resources/get_{file}.log",
	conda:
		"../envs/get_ref.yaml"
	shell:
		"""
		touch {log} && exec 2>&1 1>>{log}

		if [[ {input} == *"gz" ]]; then
			gzip -dc {input} > {output}
		else
			rsync {input} {output}
		fi
		"""


rule filter_gtf:
	input:
		txome_gtf=expand(rules.get_ref.output, file="txome.gtf", allow_missing=True),
		genome_fa=expand(rules.get_ref.output, file="genome.fa", allow_missing=True),
	output:
		gtf="{outdir}/resources/txome.filtered.gtf",
		fa="{outdir}/resources/txome.filtered.fa",
	log:
		"{outdir}/resources/filter_gtf.log",
	conda:
		"../envs/get_ref.yaml"
	shell:
		"""
		touch {log} && exec 2>&1 1>>{log}

		# String patterns used for both genomes
		ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"

		BIOTYPE_PATTERN="(protein_coding|lncRNA|lincRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
		GENE_PATTERN="gene_type \\"${{BIOTYPE_PATTERN}}\\""
		TX_PATTERN="transcript_type \\"${{BIOTYPE_PATTERN}}\\""
		READTHROUGH_PATTERN="tag \\"readthrough_transcript\\""
		PAR_PATTERN="tag \\"PAR\\""

		modified=$(mktemp)
		allowlist=$(mktemp)
		trap 'rm -f "$modified" "$allowlist"' EXIT
		cat {input.txome_gtf} | sed -E 's/gene_id "'"$ID"'";/gene_id "\\1"; gene_version "\\3";/' | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\\1"; transcript_version "\\3";/' | sed -E 's/exon_id "'"$ID"'";/exon_id "\\1"; exon_version "\\3";/' > $modified
		cat "$modified" | awk '$3 == "transcript"' | grep -E "$GENE_PATTERN" | grep -E "$TX_PATTERN" | grep -Ev "$READTHROUGH_PATTERN" | grep -Ev "$PAR_PATTERN" | sed -E 's/.*(gene_id "[^"]+").*/\\1/' | sort | uniq > $allowlist

		if [[ $(grep -E "^#" $modified | wc -l) -gt 0 ]]; then
			grep -E "^#" $modified > {output.gtf}
		fi
		grep -Ff $allowlist $modified >> {output.gtf}

		gffread {output.gtf} -g {input.genome_fa} -w {output.fa}
		"""


rule rmsk_genes:
	input:
		genome_fa=expand(rules.get_ref.output, file="genome.fa", allow_missing=True),
		rmsk_gtf=expand(rules.get_ref.output, file="rmsk.gtf", allow_missing=True),
		txome_fa=rules.filter_gtf.output.fa,
		txome_gtf=rules.filter_gtf.output.gtf,
	output:
		fa="{outdir}/resources/rmsk_genes.fa",
		gtf="{outdir}/resources/rmsk_genes.gtf",
	log:
		"{outdir}/resources/rmsk_genes.log",
	conda:
		"../envs/get_ref.yaml"
	shell:
		"""
		touch {log} && exec 2>&1 1>>{log}

		gffread -g {input.genome_fa} -w $TMP/rmsk.fa {input.rmsk_gtf}
		cat $TMP/rmsk.fa {input.txome_fa} > {output.fa}
		cat {input.rmsk_gtf} {input.txome_gtf} > {output.gtf}
		"""
