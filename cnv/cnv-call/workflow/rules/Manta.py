rule Manta_call:
	input:
		bam = bam_dir + "{sample}.ready.bam"
	output:
		vcf = TOOLS_DIR_DICT["Manta"]["output"] + "{sample}_Manta.vcf"
	params:
		vcf = TOOLS_DIR_DICT["Manta"]["info"] + "{sample}/results/variants/diploidSV.vcf.gz",
		outdir = TOOLS_DIR_DICT["Manta"]["info"] + "{sample}",
		ref = config["base"]["ref_fa"]["value"],
	conda: NGS_env_2
	threads: threads
	log: logs_dir + "Manta/{sample}.log"
	shell:
		# Configure Manta
		"(configManta.py "
		"--bam {input.bam} "
		"--referenceFasta {params.ref} "
		"--runDir {params.outdir}; "
		# Run Manta
		"{params.outdir}/runWorkflow.py --jobs {threads}; "
		"gunzip -c {params.vcf} > {output.vcf}; "
		") > {log} 2>&1"


# rule Manta_postprocess:
# 	input:
# 		vcf = rules.Manta_call.output.vcf
# 	output:
# 		vcf = TOOLS_DIR_DICT["Manta"]["output"] + "{sample}_Manta.vcf"
# 	run:
# 		Manta_postprocess(input[0], wildcards.sample, output[0])
# 		# add_var_id(output[0], "Manta")
