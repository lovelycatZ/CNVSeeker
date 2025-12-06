rule smoove_lumpy:
	input:
		bam = bam_dir + "{sample}.ready.bam"
	output:
		gz_vcf = temp(TOOLS_DIR_DICT["lumpy"]["info"] + "{sample}/{sample}-smoove.genotyped.vcf.gz"),
		vcf = protected(TOOLS_DIR_DICT["lumpy"]["output"] + "{sample}_lumpy.vcf")
	params:
		ref = REFERENCE,
		out_dir = TOOLS_DIR_DICT["lumpy"]["info"] + "{sample}/",
		exclude_regions = config["tools"]["lumpy"]["params"][f"--excl_{BUILD}"],
		exclude_chroms = "~*_hap*,~*_gl*,~_random$,~_fix$,~_alt$,~MT$,~EBV$,~^GL,hs37d5,~^chrUn,~_decoy$,~^HLA"
	conda: NGS_env_2
	benchmark: benchmark_dir + "smoove_lumpy/{sample}.log"
	log: logs_dir + "smoove_lumpy/{sample}.log"
	threads: THREADS
	shell:
		"(smoove call -x -e {params.exclude_regions} -C {params.exclude_chroms} -n {wildcards.sample} -f {params.ref} -p {threads} -o {params.out_dir} --genotype {input.bam}; "
		"gunzip -c {output.gz_vcf} > {output.vcf}; "
		") > {log} 2>&1"


# rule lumpy_postprocess:
# 	input:
# 		vcf = rules.lumpy_express.output.vcf
# 	output:
# 		vcf = TOOLS_DIR_DICT["lumpy"]["output"] + "{sample}_lumpy.vcf"
# 	run:
# 		lumpy_postprocess(input[0], wildcards.sample, output[0])
# 		# add_var_id(output[0], "lumpy")
