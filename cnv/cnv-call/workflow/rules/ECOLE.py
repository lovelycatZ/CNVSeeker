rule create_read_depth:
	input:
		bam = bam_dir + "{sample}.ready.bam"
	output:
		rd = TOOLS_DIR_DICT["ECOLE"]["info"] + "read_depth_{sample}/{sample}.rd.txt"
	params:
		bed = EXOME_INTERVALS
	conda: NGS_env_1
	benchmark: benchmark_dir + "ECOLE/read_depth/{sample}.log"
	threads: THREADS
	log: logs_dir + "ECOLE/read_depth/{sample}.log"
	shell:
		"(sambamba depth base -L {params.bed} {input.bam} > {output.rd}; "
		") > {log} 2>&1"


rule ECOLE_preprocess_sample:
	input:
		rd = TOOLS_DIR_DICT["ECOLE"]["info"] + "read_depth_{sample}/{sample}.rd.txt"
	output:
		labeled_data = TOOLS_DIR_DICT["ECOLE"]["info"] + "processed_sample_{sample}/{sample}_labeled_data.npy"
	params:
		script = config["tools"]["ECOLE"]["bin"]["--preprocess_sample"],
		bed = EXOME_INTERVALS,
	conda: ECOLE
	benchmark: benchmark_dir + "ECOLE/preprocess_sample/{sample}.log"
	threads: THREADS
	log: logs_dir + "ECOLE/preprocess_sample/{sample}.log"
	shell:
		"(python {params.script} -rd {input.rd} -t {params.bed} -o {output.labeled_data}; "
		") > {log} 2>&1"


rule ECOLE_call:
	input:
		labeled_data = TOOLS_DIR_DICT["ECOLE"]["info"] + "processed_sample_{sample}/{sample}_labeled_data.npy"
	output:
		call = TOOLS_DIR_DICT["ECOLE"]["info"] + "call_{sample}/{sample}.csv"
	params:
		script = config["tools"]["ECOLE"]["bin"]["--ecole_call"],
		bed = EXOME_INTERVALS,
		normalize_data = config["tools"]["ECOLE"]["params"]["--normalize_data"],
		labeled_data_dir = TOOLS_DIR_DICT["ECOLE"]["info"] + "processed_sample_{sample}",
		call_dir = TOOLS_DIR_DICT["ECOLE"]["info"] + "call_{sample}",
		tmp1 = TOOLS_DIR_DICT["ECOLE"]["info"] + "tmp1_{sample}",
		tmp2 = TOOLS_DIR_DICT["ECOLE"]["info"] + "tmp2_{sample}",
	conda: ECOLE
	benchmark: benchmark_dir + "ECOLE/call/{sample}.log"
	threads: THREADS
	log: logs_dir + "ECOLE/call/{sample}.log"
	shell:
		"(python {params.script} --model ecole --input {params.labeled_data_dir} --output {params.call_dir} --cnv merged --batch_size 10 --normalize {params.normalize_data} --gpu {threads} -t1 {params.tmp1} -t2 {params.tmp2} -t {params.bed}; "
		") > {log} 2>&1"


rule ECOLE_convert:
	input:
		call = TOOLS_DIR_DICT["ECOLE"]["info"] + "call_{sample}/{sample}.csv"
	output:
		vcf = TOOLS_DIR_DICT["ECOLE"]["output"] + "{sample}_ECOLE.vcf"
	params:
		script = HOME_DIR + "workflow/scripts/ecole2vcf.py"
	conda: NGS_env_1
	benchmark: benchmark_dir + "ECOLE/convert/{sample}.log"
	threads: THREADS
	log: logs_dir + "ECOLE/convert/{sample}.log"
	shell:
		"(python {params.script} {input.call} {wildcards.sample} {output.vcf}; "
		") > {log} 2>&1"
