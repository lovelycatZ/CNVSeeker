tool = os.path.basename(workflow.snakefile).split(".py")[0]
cuteSV_output_dir, cuteSV_info_dir = TOOLS_DIR_DICT[tool].values()



rule cuteSV:
	input:
		# bam = lambda wildcards: bam_dir + "{sample}.{platform}.{aligner}.ready.bam".format(sample = wildcards.sample, platform = get_platform(wildcards), aligner = "minimap2")
		bam = lambda wildcards: get_caller_input_sample(wildcards, os.path.basename(workflow.snakefile).split(".py")[0])
	output:
		vcf = cuteSV_output_dir + "{{sample}}_{tool}.vcf".format(tool = os.path.basename(workflow.snakefile).split(".py")[0])
	params:
		pt_specific_params = lambda wildcards: get_cuteSV_params(wildcards),
		caller_filter_params = get_caller_filter_params(os.path.basename(workflow.snakefile).split(".py")[0], config["tools"]["TGS_caller_filter_options"]["params"]),
		ref = REFERENCE,
		work_dir = cuteSV_info_dir + "{sample}",
		sample_name = lambda wildcards: "{sample}_{tool}".format(sample = wildcards.sample, tool = os.path.basename(workflow.snakefile).split(".py")[0])
	conda: TGS_env_1
	benchmark: benchmark_dir + "cuteSV/{sample}.log"
	threads: THREADS
	log: logs_dir + "cuteSV/{sample}.log"
	shell:
		"(mkdir -p {params.work_dir}; "
		"cuteSV {params.pt_specific_params} {params.caller_filter_params} "
		"-t {threads} -S {wildcards.sample} --genotype {input.bam} {params.ref} {output.vcf} {params.work_dir}; "
		"sed -i '/^#CHROM/s/\t[^\t]*$/\t{params.sample_name}/' {output.vcf}; "
		") > {log} 2>&1"
