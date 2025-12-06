tool = os.path.basename(workflow.snakefile).split(".py")[0]
SVDSS_output_dir, SVDSS_info_dir = TOOLS_DIR_DICT[tool].values()



rule SVDSS_index:
	input:
		ref = REFERENCE
	output:
		fmd = SVDSS_info_dir + "ref.fmd"
	params:
		SVDSS = config["tools"]["SVDSS"]["bin"]
	conda: TGS_env_1
	benchmark: benchmark_dir + "SVDSS/index.log"
	threads: THREADS
	log: logs_dir + "SVDSS/index.log"
	wildcard_constraints:
		sample = "|".join(fa_samples) if fa_samples else "no_wildcards"	
	shell:
		"({params.SVDSS} index --reference {input.ref} --index {output.fmd}; "
		") > {log} 2>&1"


rule SVDSS_smooth:
	input:
		bam = lambda wildcards: get_caller_input_sample(wildcards, os.path.basename(workflow.snakefile).split(".py")[0])
	output:
		bam = SVDSS_info_dir + "{sample}/{sample}.smoothed.bam",
		bai = SVDSS_info_dir + "{sample}/{sample}.smoothed.bam.bai"
	params:
		ref = REFERENCE
		SVDSS = config["tools"]["SVDSS"]["bin"]
	conda: TGS_env_1
	benchmark: benchmark_dir + "SVDSS/smooth/{sample}.log"
	threads: THREADS
	log: logs_dir + "SVDSS/smooth/{sample}.log"
	wildcard_constraints:
		sample = "|".join(fa_samples) if fa_samples else "no_wildcards"	
	shell:
		"({params.SVDSS} smooth --reference {params.ref} --bam {input.bam} --threads {threads} > {output.bam}; "
		"samtools index -@ {threads} {output.bam}; "
		") > {log} 2>&1"


rule SVDSS_search:
	input:
		bam = SVDSS_info_dir + "{sample}/{sample}.smoothed.bam",
		fmd = SVDSS_info_dir + "ref.fmd"
	output:
		sfs = SVDSS_info_dir + "{sample}/specifics.txt"
	params:
		SVDSS = config["tools"]["SVDSS"]["bin"]
	conda: TGS_env_1
	benchmark: benchmark_dir + "SVDSS/search/{sample}.log"
	threads: THREADS
	log: logs_dir + "SVDSS/search/{sample}.log"
	wildcard_constraints:
		sample = "|".join(fa_samples) if fa_samples else "no_wildcards"	
	shell:
		"({params.SVDSS} search --index {input.fmd} --bam {input.bam} --threads {threads} > {output.sfs}; "
		") > {log} 2>&1"


rule SVDSS_call:
	input:
		bam = SVDSS_info_dir + "{sample}/{sample}.smoothed.bam",
		sfs = SVDSS_info_dir + "{sample}/specifics.txt"
	output:
		vcf = SVDSS_output_dir + "{{sample}}_{tool}.vcf".format(tool = os.path.basename(workflow.snakefile).split(".py")[0])
	params:
		caller_filter_params = get_caller_filter_params(os.path.basename(workflow.snakefile).split(".py")[0], config["tools"]["TGS_caller_filter_options"]["params"]),
		sample_name = lambda wildcards: "{sample}_{tool}".format(sample = wildcards.sample, tool = os.path.basename(workflow.snakefile).split(".py")[0]),
		ref = REFERENCE,
		SVDSS = config["tools"]["SVDSS"]["bin"]
	conda: TGS_env_1
	benchmark: benchmark_dir + "SVDSS/call/{sample}.log"
	threads: THREADS
	log: logs_dir + "SVDSS/call/{sample}.log"
	wildcard_constraints:
		sample = "|".join(fa_samples) if fa_samples else "no_wildcards"
	shell:
		"({params.SVDSS} call {params.caller_filter_params} --reference {params.ref} --bam {input.bam} "
		"--sfs {input.sfs} --threads {threads} | bcftools sort > {output.vcf}; "
		"sed -i '/^#CHROM/s/\t[^\t]*$/\t{params.sample_name}/' {output.vcf}; "
		") > {log} 2>&1"
