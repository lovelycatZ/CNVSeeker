rule CNVpytor_call:
	input:
		bam = bam_dir + "{sample}.ready.bam"
	output:
		pytor = temp(TOOLS_DIR_DICT["CNVpytor"]["info"] + "{sample}/{sample}.pytor"),
		tsv = TOOLS_DIR_DICT["CNVpytor"]["info"] + "{sample}/{sample}.tsv",
		# vcf   = TOOLS_DIR_DICT["CNVpytor"]["info"] + "{sample}/{sample}.vcf"
		vcf = protected(TOOLS_DIR_DICT["CNVpytor"]["output"] + "{sample}_CNVpytor.vcf")
	params:
		script		= config["tools"]["CNVpytor"]["params"]["--convertor"],
		window_size = config["tools"]["CNVpytor"]["params"]["--window_size"],
		chroms 		= " ".join(CHROM_LIST),
		ref 		= REFERENCE
	conda: NGS_env_1
	threads: THREADS
	benchmark: benchmark_dir + "CNVpytor/call/{sample}.log"
	log: logs_dir + "CNVpytor/call/{sample}.log"
	shell:
		"(cnvpytor -root {output.pytor} -rd {input.bam} -chrom {params.chroms}; "
		"cnvpytor -root {output.pytor} -his {params.window_size} -chrom {params.chroms}; "
		"cnvpytor -root {output.pytor} -partition {params.window_size}; "
		"cnvpytor -root {output.pytor} -call {params.window_size} > {output.tsv}; "
		"python {params.script} {output.tsv} {wildcards.sample} {output.vcf}; "
		") > {log} 2>&1"


# rule CNVpytor_postprocess:
# 	input:
# 		vcf  = rules.CNVpytor_call.output.vcf
# 	output:
# 		vcf = TOOLS_DIR_DICT["CNVpytor"]["output"] + "{sample}_CNVpytor.vcf"
# 	params:
# 		p1_thd = str(config["tools"]["CNVpytor"]["params"]["--p1"]),
# 		q0_thd = str(config["tools"]["CNVpytor"]["params"]["--q0"])
# 	run:
# 		CNVpytor_postprocess(input[0], wildcards.sample, output[0], params["p1_thd"], params["q0_thd"])
# 		# add_var_id(output[0], "CNVpytor")
