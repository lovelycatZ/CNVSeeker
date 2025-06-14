rule cnmops_call:
	input:
		bams = expand(bam_dir + "{sample}.ready.bam", sample = SAMPLES)
	output:
		cnvs = TOOLS_DIR_DICT["cn.mops"]["info"] + "cnvs.tsv",
		cnvr = TOOLS_DIR_DICT["cn.mops"]["info"] + "cnvr.tsv"
	params:
		script 			= HOME_DIR + "workflow/scripts/cn.mops.R",
		bed 			= EXOME_INTERVALS,
		window_size 	= config["tools"]["cn_mops"]["params"]["--window_size"],
		chrom_style 	= CHROM_STYLE,
		ngs_data_type 	= NGS_DATA_TYPE
	conda: NGS_env_1
	benchmark: benchmark_dir + f"cn.mops/call_{NGS_DATA_TYPE}.log"
	log: logs_dir + f"cn.mops/call_{NGS_DATA_TYPE}.log"
	threads: THREADS
	shell:
		"(Rscript {params.script} {params.ngs_data_type} {params.chrom_style} {params.window_size} {params.bed} {threads} {output.cnvs} {output.cnvr} {input.bams}; "
		") > {log} 2>&1"


rule cnmops_postprocess:
	input:
		cnvs = TOOLS_DIR_DICT["cn.mops"]["info"] + "cnvs.tsv"
	output:
		vcf = protected(TOOLS_DIR_DICT["cn.mops"]["output"] + "{sample}_cn.mops.vcf")
	benchmark: benchmark_dir + "cn.mops/postprocess/{sample}.log"
	run:
		cnmops_postprocess(input[0], wildcards.sample, output[0])
		# add_var_id(output[0], "cn.mops")
