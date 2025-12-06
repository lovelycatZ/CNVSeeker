import os


REFERENCE_DIR = os.path.dirname(REFERENCE) + "/"

rule ControlFREEC_chr_files:
	input:
		ref = REFERENCE
	output:
		chr_file = REFERENCE_DIR + "{chr}.fa"
	params:
		chr = "{chr}"
	conda: NGS_env_1
	threads: THREADS
	benchmark: benchmark_dir + "ControlFREEC/chr_files/{chr}.log"
	log: logs_dir + "ControlFREEC/chr_files/{chr}.log"
	shell:
		"(samtools faidx {input.ref} {params.chr} > {output.chr_file}; "
		") > {log} 2>&1"



rule ControlFREEC_new_fai:
	input:
		chrlen = config["tools"]["ControlFREEC"]["params"]["--chr_len"]
	output:
		chrlen = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "chrlen.new.txt"
	params:
		chrom_list = CHROM_LIST
	run:
		with open(input["chrlen"], "r") as f, open(output["chrlen"], "w") as out:
			for line in f:
				chr = line.strip().split("\t")[0]
				if chr not in params["chrom_list"]:
					continue
				out.write(line)



rule ControlFREEC_config:
	input:
		chrlen = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "chrlen.new.txt",
		chr_files = expand(REFERENCE_DIR + "{chr}.fa", chr = CHROM_LIST),
		bam_sample = bam_dir + "{sample}.ready.bam",
	output:
		config = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "config/{sample}.config.txt"
	params:
		script = HOME_DIR + "workflow/scripts/freec_config.py",
		bed = EXOME_INTERVALS,
		chrdir = REFERENCE_DIR,
		chrlen = config["tools"]["ControlFREEC"]["params"]["--chr_len"],
		window_size = config["tools"]["ControlFREEC"]["params"]["--window_size"],
		sambamba = config["envs"]["conda"]["conda_env_dir"] + f"/{NGS_env_1}/bin/sambamba",
		out_dir = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "temp/",
		sex = lambda wildcards: get_sex(wildcards)
	conda: NGS_env_1
	threads: THREADS
	benchmark: benchmark_dir + "ControlFREEC/make_config/{sample}.log"
	log: logs_dir + "ControlFREEC/make_config/{sample}.log"
	shell:
		"(python {params.script} "
		# general
		"--outputDir {params.out_dir} "
		"--chrFiles {params.chrdir} "
		"--chrLenFile {input.chrlen} "
		"--window {params.window_size} "
		"--ploidy 2 "
		"--sex {params.sex} "
		"--breakPointType 2 "
		"--breakPointThreshold 0.8 "
		"--maxThreads {threads} "
		"--sambamba {params.sambamba} "
		"--SambambaThreads {threads} "
		"--noisyData TRUE "
		"--printNA FALSE "
		# sample
		"--input_sample {input.bam_sample} "
		"--input_sample_format BAM "
		"--input_sample_orientation FR "
		# end
		"--out_config {output.config}; "
		") > {log} 2>&1"
		

rule ControlFREEC_call:
	input:
		config = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "config/{sample}.config.txt"
	output:
		cnv   = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "temp/{sample}.ready.bam_CNVs",
		ratio = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "temp/{sample}.ready.bam_ratio.txt"
	conda: NGS_env_1
	threads: THREADS
	benchmark: benchmark_dir + "ControlFREEC/call/{sample}.log"
	log: logs_dir + "ControlFREEC/call/{sample}.log"
	shell:
		"(freec -conf {input.config}; "
		") > {log} 2>&1"


rule ControlFREEC_significance:
	input:
		cnv = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "temp/{sample}.ready.bam_CNVs",
		ratio = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "temp/{sample}.ready.bam_ratio.txt"
	output:
		pvalue_txt = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "temp/{sample}.ready.bam_CNVs.p.value.txt"
	params:
		script = HOME_DIR + "workflow/scripts/assess_significance.R"
	conda: NGS_env_1
	threads: THREADS
	benchmark: benchmark_dir + "ControlFREEC/significance/{sample}.log"
	log: logs_dir + "ControlFREEC/significance/{sample}.log"
	shell:
		"(Rscript {params.script} {input.cnv} {input.ratio}; "
		") > {log} 2>&1"


rule ControlFREEC_convert:
	input:
		pvalue_txt = TOOLS_DIR_DICT["ControlFREEC"]["info"] + "temp/{sample}.ready.bam_CNVs.p.value.txt"
	output:
		vcf = protected(TOOLS_DIR_DICT["ControlFREEC"]["output"] + "{sample}_ControlFREEC.vcf")
	params:
		script = HOME_DIR + "workflow/scripts/freec2vcf.py",
		chrom_style = CHROM_STYLE
	conda: NGS_env_1
	threads: THREADS
	benchmark: benchmark_dir + "ControlFREEC/convert/{sample}.log"
	log: logs_dir + "ControlFREEC/convert/{sample}.log"
	shell:
		"(python {params.script} {input.pvalue_txt} {wildcards.sample} {params.chrom_style} {output.vcf}; "
		") > {log} 2>&1"
