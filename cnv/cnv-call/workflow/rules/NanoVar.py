tool = os.path.basename(workflow.snakefile).split(".py")[0]
NanoVar_output_dir, NanoVar_info_dir = TOOLS_DIR_DICT[tool].values()



rule NanoVar:
	input:
		# bam = lambda wildcards: bam_dir + "{sample}.{platform}.{aligner}.ready.bam".format(sample = wildcards.sample, platform = get_platform(wildcards), aligner = "minimap2")
		bam = lambda wildcards: get_caller_input_sample(wildcards, os.path.basename(workflow.snakefile).split(".py")[0])
	output:
		vcf  = NanoVar_output_dir + "{{sample}}_{tool}.vcf".format(tool = os.path.basename(workflow.snakefile).split(".py")[0])
	params:
		caller_filter_params = get_caller_filter_params(os.path.basename(workflow.snakefile).split(".py")[0], config["tools"]["TGS_caller_filter_options"]["params"]),
		vcf = lambda wildcards: get_Nanovar_out_vcf(wildcards, NanoVar_info_dir),
		pt_specific_params = lambda wildcards: get_NanoVar_params(wildcards),
		ref = REFERENCE,
		work_dir = NanoVar_info_dir + "{sample}",
		sample_name = lambda wildcards: "{sample}_{tool}".format(sample = wildcards.sample, tool = os.path.basename(workflow.snakefile).split(".py")[0])
	conda: TGS_env_1
	benchmark: benchmark_dir + "NanoVar/{sample}.log"
	threads: THREADS
	log: logs_dir + "NanoVar/{sample}.log"
	shell:
		"(nanovar {params.pt_specific_params} {params.caller_filter_params} "
		"-t {threads} -f hg38 {input.bam} {params.ref} {params.work_dir}; "
		"cp {params.vcf} {output.vcf}; "
		"sed -i '/^#CHROM/s/\t[^\t]*$/\t{params.sample_name}/' {output.vcf}; "
		") > {log} 2>&1"
