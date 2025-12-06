REF_FAI = REFERENCE + ".fai"

if not os.path.exists(REF_FAI):
	rule gridss_reference:
		input:
			ref = REFERENCE
		output:
			ref_fai = REF_FAI
		params:
			jar = config["bin"]["gridss_jar"]
		conda: NGS_env_1
		benchmark: benchmark_dir + "gridss/reference.log"
		log: logs_dir + "gridss/reference.log"
		threads: THREADS
		shell:
			"(gridss --jvmheap 15g --otherjvmheap 15g -r {input.ref} -j {params.jar} -t {threads} -s setupreference; "
			") > {log} 2>&1"


rule gridss_preprocess:
	input:
		bam = bam_dir + "{sample}.ready.bam",
		ref_fai = REF_FAI,
	output:
		sv_bam = temp(TOOLS_DIR_DICT["gridss"]["info"] + "{sample}/{sample}.ready.bam.gridss.working/{sample}.ready.bam.sv.bam")
	params:
		jar = config["bin"]["gridss_jar"],
		ref = REFERENCE,
		work_dir = TOOLS_DIR_DICT["gridss"]["info"] + "{sample}/"
	conda: NGS_env_1
	benchmark: benchmark_dir + "gridss/preprocess/{sample}.log"
	log: logs_dir + "gridss/preprocess/{sample}.log"
	threads: THREADS
	shell:
		"(cd {params.work_dir}; "
		"gridss --jvmheap 15g --otherjvmheap 15g -r {params.ref} -j {params.jar} -t {threads} -s preprocess {input.bam}; "
		") > {log} 2>&1"


rule gridss_assemble:
	input:
		bam = bam_dir + "{sample}.ready.bam",
		sv_bam = TOOLS_DIR_DICT["gridss"]["info"] + "{sample}/{sample}.ready.bam.gridss.working/{sample}.ready.bam.sv.bam"
	output:
		assemble_bam = temp(TOOLS_DIR_DICT["gridss"]["info"] + "{sample}.assembly.bam")
	params:
		jar = config["bin"]["gridss_jar"],
		ref = REFERENCE,
		work_dir = TOOLS_DIR_DICT["gridss"]["info"] + "{sample}/"
	conda: NGS_env_1
	benchmark: benchmark_dir + "gridss/assemble/{sample}.log"
	log: logs_dir + "gridss/assemble/{sample}.log"
	threads: THREADS
	shell:
		"(cd {params.work_dir}; "
		"gridss --jvmheap 15g --otherjvmheap 15g -r {params.ref} -j {params.jar} -t {threads} -s assemble -a {output.assemble_bam} {input.bam}; "
		") > {log} 2>&1"


rule gridss_call:
	input:
		bam = bam_dir + "{sample}.ready.bam",
		sv_bam = TOOLS_DIR_DICT["gridss"]["info"] + "{sample}/{sample}.ready.bam.gridss.working/{sample}.ready.bam.sv.bam",
		assemble_bam = TOOLS_DIR_DICT["gridss"]["info"] + "{sample}.assembly.bam"
	output:
		vcf = TOOLS_DIR_DICT["gridss"]["output"] + "{sample}_gridss.vcf"
	params:
		jar = config["bin"]["gridss_jar"],
		ref = REFERENCE,
		work_dir = TOOLS_DIR_DICT["gridss"]["info"] + "{sample}/"
	conda: NGS_env_1
	benchmark: benchmark_dir + "gridss/call/{sample}.log"
	log: logs_dir + "gridss/call/{sample}.log"
	threads: THREADS
	shell:
		"(cd {params.work_dir}; "
		"gridss --jvmheap 15g --otherjvmheap 15g -r {params.ref} -j {params.jar} -t {threads} -s call -a {input.assemble_bam} -o {output.vcf} {input.bam}; "
		") > {log} 2>&1"
