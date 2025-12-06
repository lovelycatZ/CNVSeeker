rule delly_call:
	input:
		bam = bam_dir + "{sample}.ready.bam"
	output:
		vcf = protected(TOOLS_DIR_DICT["delly"]["output"] + "{sample}_delly.vcf")
	params:
		ref = REFERENCE,
		excl = config["tools"]["delly"]["params"][f"--excl_{BUILD}"]
	conda: NGS_env_1
	benchmark: benchmark_dir + "delly/{sample}.log"
	log: logs_dir + "delly/{sample}.log"
	threads: THREADS
	shell:
		"(delly call -x {params.excl} -g {params.ref} -t DUP,DEL -z 4 -s 12 {input.bam} > {output.vcf}; "
		") > {log} 2>&1"


# rule delly_postprocess:
# 	input:
# 		vcf = rules.concat_vcf.output.vcf
# 	output:
# 		vcf = TOOLS_DIR_DICT["delly"]["output"] + "{sample}_delly.vcf"
# 	run:
# 		delly_postprocess(input[0], wildcards.sample, output[0])
# 		# add_var_id(output[0], "delly")
