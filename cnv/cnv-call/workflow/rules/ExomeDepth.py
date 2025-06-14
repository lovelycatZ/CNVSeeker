rule ExomeDepth_call:
	input:
		bams = expand(bam_dir + "{sample}.ready.bam", sample = SAMPLES),
	output:
		tsv = expand(TOOLS_DIR_DICT["ExomeDepth"]["info"] + "{sample}.tsv", sample = SAMPLES)
	params:
		infos_dir 	= TOOLS_DIR_DICT["ExomeDepth"]["info"],
		exomedepth 	= HOME_DIR + "/workflow/scripts/ExomeDepth.R",
		ref 		= REFERENCE,
		bed 		= EXOME_INTERVALS
	conda: NGS_env_1
	benchmark: benchmark_dir + "ExomeDepth/call.log"
	threads: THREADS
	log: logs_dir + "ExomeDepth/call.log"
	shell:
		"(Rscript {params.exomedepth} {params.ref} {params.bed} {params.infos_dir} {input.bams}; "
		") > {log} 2>&1"


rule ExomeDepth_postprocess:
	input:
		tsv = TOOLS_DIR_DICT["ExomeDepth"]["info"] + "{sample}.tsv"
	output:
		vcf = TOOLS_DIR_DICT["ExomeDepth"]["output"] + "{sample}_ExomeDepth.vcf"
	benchmark: benchmark_dir + "ExomeDepth/postprocess/{sample}.log"
	run:
		ExomeDepth_postprocess(input[0], wildcards.sample, output[0])
		# add_var_id(output[0], "ExomeDepth")
