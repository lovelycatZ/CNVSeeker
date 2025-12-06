tool = os.path.basename(workflow.snakefile).split(".py")[0]
pbsv_output_dir, pbsv_info_dir = TOOLS_DIR_DICT[tool].values()



rule pbsv_discovery:
	input:
		# bam = lambda wildcards: bam_dir + "{sample}.{platform}.{aligner}.ready.bam".format(sample = wildcards.sample, platform = get_platform(wildcards), aligner = "pbmm2")
		bam = lambda wildcards: get_caller_input_sample(wildcards, os.path.basename(workflow.snakefile).split(".py")[0])
	output:
		svsig = pbsv_info_dir + "{sample}.svsig.gz"
	params:
		ref = REFERENCE
	conda: TGS_env_1
	benchmark: benchmark_dir + "pbsv/discovery/{sample}.log"
	threads: THREADS
	log: logs_dir + "pbsv/discovery/{sample}.log"
	shell:
		"(pbsv discover {input.bam} {output.svsig}; "
		") > {log} 2>&1"


rule pbsv_call:
	input:
		svsig = pbsv_info_dir + "{sample}.svsig.gz"
	output:
		vcf = pbsv_output_dir + "{{sample}}_{tool}.vcf".format(tool = os.path.basename(workflow.snakefile).split(".py")[0])
	params:
		caller_filter_params = get_caller_filter_params(os.path.basename(workflow.snakefile).split(".py")[0], config["tools"]["TGS_caller_filter_options"]["params"]),
		ref = REFERENCE,
		sample_name = lambda wildcards: "{sample}_{tool}".format(sample = wildcards.sample, tool = os.path.basename(workflow.snakefile).split(".py")[0])
	conda: TGS_env_1
	benchmark: benchmark_dir + "pbsv/call/{sample}.log"
	threads: THREADS
	log: logs_dir + "pbsv/call/{sample}.log"
	shell:
		"(pbsv call {params.caller_filter_params} --num-threads {threads} {params.ref} {input.svsig} {output.vcf}; "
		"sed -i '/^#CHROM/s/\t[^\t]*$/\t{params.sample_name}/' {output.vcf}; "
		") > {log} 2>&1"
