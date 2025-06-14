tool = os.path.basename(workflow.snakefile).split(".py")[0]
SVision_pro_output_dir, SVision_pro_info_dir = TOOLS_DIR_DICT[tool].values()



rule SVision_pro:
	input:
		# bam = lambda wildcards: bam_dir + "{sample}.{platform}.{aligner}.ready.bam".format(sample = wildcards.sample, platform = get_platform(wildcards), aligner = "winnowmap")
		bam = lambda wildcards: get_caller_input_sample(wildcards, os.path.basename(workflow.snakefile).split(".py")[0])
	output:
		vcf  = SVision_pro_output_dir + "{{sample}}_{tool}.vcf".format(tool = os.path.basename(workflow.snakefile).split(".py")[0]),
		sample_dir = directory(SVision_pro_info_dir + "{sample}")
	params:
		caller_filter_params = get_caller_filter_params(os.path.basename(workflow.snakefile).split(".py")[0], config["tools"]["TGS_caller_filter_options"]["params"]),
		vcf = lambda wildcards: SVision_pro_info_dir + "{sample}/{sample}.svision_pro_v1.8.s{min_support}.vcf".format(sample = wildcards.sample, min_support = config["tools"]["TGS_caller_filter_options"]["params"]["--min_support"]),
		access_bed = config["tools"]["SVision_pro"]["params"]["--access_bed"],
		model = config["tools"]["SVision_pro"]["params"]["--model"],
		bin = config["tools"]["SVision_pro"]["bin"],
		ref = REFERENCE,
		sample_name = lambda wildcards: "{sample}_{tool}".format(sample = wildcards.sample, tool = os.path.basename(workflow.snakefile).split(".py")[0])
	conda: TGS_env_2
	benchmark: benchmark_dir + "SVision_pro/{sample}.log"
	threads: THREADS
	log: logs_dir + "SVision_pro/{sample}.log"
	shell:
		"({params.bin} {params.caller_filter_params} "
		"--target_path {input.bam} "
		"--genome_path {params.ref} "
		"--model_path {params.model} "
		"--access_path {params.access_bed} "
		"--sample {wildcards.sample} "
		"--out_path {output.sample_dir} "
		"--detect_mode germline "
		"--process {threads}; "
		"cp {params.vcf} {output.vcf}; "
		"sed -i '/^#CHROM/s/\t[^\t]*$/\t{params.sample_name}/' {output.vcf}; "
		") > {log} 2>&1"