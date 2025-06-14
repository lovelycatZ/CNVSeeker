rule Wham_call:
	input:
		bam = bam_dir + "{sample}.ready.bam"
	output:
		vcf = protected(TOOLS_DIR_DICT["Wham"]["output"] + "{sample}_Wham.vcf")
	params:
		ref = REFERENCE,
		chrom_list = CHROM_LIST
	conda: NGS_env_1
	benchmark: benchmark_dir + "Wham/{sample}.log"
	log: logs_dir + "Wham/{sample}.log"
	threads: THREADS
	shell:
		"(whamg -f {input.bam} -a {params.ref} -c {params.chrom_list} -x {threads} > {output.vcf}; "
		") > {log} 2>&1"
