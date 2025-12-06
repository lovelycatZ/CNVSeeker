tool = os.path.basename(workflow.snakefile).split(".py")[0]
DeBreak_output_dir, DeBreak_info_dir = TOOLS_DIR_DICT[tool].values()



rule DeBreak:
	input:
		# bam = lambda wildcards: bam_dir + "{sample}.{platform}.{aligner}.ready.bam".format(sample = wildcards.sample, platform = get_platform(wildcards), aligner = "minimap2")
		bam = lambda wildcards: get_caller_input_sample(wildcards, os.path.basename(workflow.snakefile).split(".py")[0])
	output:
		vcf  = DeBreak_output_dir + "{{sample}}_{tool}.vcf".format(tool = os.path.basename(workflow.snakefile).split(".py")[0]),
		sample_dir = directory(DeBreak_info_dir + "{sample}")
	params:
		caller_filter_params = get_caller_filter_params(os.path.basename(workflow.snakefile).split(".py")[0], config["tools"]["TGS_caller_filter_options"]["params"]),
		vcf = DeBreak_info_dir + "{sample}/debreak.vcf",
		ref = REFERENCE,
		sample_name = lambda wildcards: "{sample}_{tool}".format(sample = wildcards.sample, tool = os.path.basename(workflow.snakefile).split(".py")[0])
	conda: TGS_env_2
	benchmark: benchmark_dir + "DeBreak/{sample}.log"
	threads: THREADS
	log: logs_dir + "DeBreak/{sample}.log"
	shell:
		"(debreak {params.caller_filter_params} -t {threads} --rescue_dup --poa "
		"--bam {input.bam} -o {output.sample_dir} -r {params.ref}; "
		"cp {params.vcf} {output.vcf}; "
		"sed -i '/^#CHROM/s/\t[^\t]*$/\t{params.sample_name}/' {output.vcf}; "
		") > {log} 2>&1"