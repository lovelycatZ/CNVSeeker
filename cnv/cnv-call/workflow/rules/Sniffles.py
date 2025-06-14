tool = os.path.basename(workflow.snakefile).split(".py")[0]
Sniffles_output_dir, Sniffles_info_dir = TOOLS_DIR_DICT[tool].values()



rule Sniffles:
	input:
		# bam = lambda wildcards: bam_dir + "{sample}.{platform}.{aligner}.ready.bam".format(sample = wildcards.sample, platform = get_platform(wildcards), aligner = "winnowmap")
		bam = lambda wildcards: get_caller_input_sample(wildcards, os.path.basename(workflow.snakefile).split(".py")[0])
	output:
		vcf  = Sniffles_output_dir + "{{sample}}_{tool}.vcf".format(tool = os.path.basename(workflow.snakefile).split(".py")[0])
	params:
		caller_filter_params = get_caller_filter_params(os.path.basename(workflow.snakefile).split(".py")[0], config["tools"]["TGS_caller_filter_options"]["params"]),
		tandem_repeats = config["tools"]["Sniffles"]["params"][f"--tandem_repeats_{BUILD}"],
		ref = REFERENCE,
		sample_name = lambda wildcards: "{sample}_{tool}".format(sample = wildcards.sample, tool = os.path.basename(workflow.snakefile).split(".py")[0])
	conda: TGS_env_1
	benchmark: benchmark_dir + "Sniffles/{sample}.log"
	threads: THREADS
	log: logs_dir + "Sniffles/{sample}.log"
	shell:
		"(sniffles {params.caller_filter_params} --tandem-repeats {params.tandem_repeats} --threads {threads} "
		"--input {input.bam} --vcf {output.vcf} --reference {params.ref}; "
		"sed -i '/^#CHROM/s/\t[^\t]*$/\t{params.sample_name}/' {output.vcf}; "
		") > {log} 2>&1"
