import os


def get_final_out():
	out_list = []
	if INTEGRATE:
		out_dir = merged_vcf_dir + "output/"
	else:
		out_dir = final_dir
	
	for sample in SAMPLES_S:
		out_file = out_dir + f"{sample}.{CALL_OUTPUT_FORMAT}"
		out_list.append(out_file)
		
	return out_list


def get_fastq_sample(wildcards, orientation):
	sample = wildcards.sample
	if DF_FQ.loc[sample, 'fq_r1'].endswith((".fq", ".fastq")):
		return fastq_raw_dir + f"{sample}_{orientation}.fq"
	else:
		return fastq_raw_dir + f"{sample}_{orientation}.fq.gz"
	

def create_read_count_str(file_dir, suffix):
	lst = []
	for sample in SAMPLES:
		file = file_dir + f"{sample}.{suffix}"
		lst.append(f"-I {file}")
	return " ".join(lst)


def create_calls_shard_str(file_dir, indexs):
	lst = []
	for index in indexs:
		file = file_dir + f"cohort_{index}-calls"
		lst.append(f"--calls-shard-path {file}")
	return " ".join(lst)


def create_model_shard_str(file_dir, indexs):
	lst = []
	for index in indexs:
		file = file_dir + f"cohort_{index}-model"
		lst.append(f"--model-shard-path {file}")
	return " ".join(lst)


def get_sex(wildcards):
	if wildcards.sample in DF_FQ.index:
		sex = DF_FQ.loc[wildcards.sample, 'sex']
	else:
		sex = DF_BAM.loc[wildcards.sample, 'sex']
	
	if sex == "F":
		return	"XX"
	else:
		return "XY"
	

rule link_fastq:
	input:
		fastq_r1 = lambda wildcards: DF_FQ.loc[DF_FQ.index == wildcards.sample, 'fq_r1'],
		fastq_r2 = lambda wildcards: DF_FQ.loc[DF_FQ.index == wildcards.sample, 'fq_r2']
	output:
		ln_fastq_r1 = fastq_raw_dir + "{sample}_R1.fq",
		ln_fastq_r2 = fastq_raw_dir + "{sample}_R2.fq"
	wildcard_constraints:
		sample = "|".join(fq_samples) if fq_samples else "no_wildcards"
	run:
		for fq, ln_fq in zip(input, output):
			if not os.path.islink(ln_fq):
				os.symlink(fq, ln_fq)


rule link_fastq_gz:
	input:
		fastq_r1 = lambda wildcards: DF_FQ.loc[DF_FQ.index == wildcards.sample, 'fq_r1'],
		fastq_r2 = lambda wildcards: DF_FQ.loc[DF_FQ.index == wildcards.sample, 'fq_r2']
	output:
		ln_fastq_r1 = fastq_raw_dir + "{sample}_R1.fq.gz",
		ln_fastq_r2 = fastq_raw_dir + "{sample}_R2.fq.gz"
	wildcard_constraints:
		sample = "|".join(fq_gz_samples) if fq_gz_samples else "no_wildcards"
	run:
		for fq, ln_fq in zip(input, output):
			if not os.path.islink(ln_fq):
				os.symlink(fq, ln_fq)

				
rule link_bam:
	input:
		bam = lambda wildcards: DF_BAM.loc[DF_BAM.index == wildcards.sample, 'bam'],
		bai = lambda wildcards: DF_BAM.loc[DF_BAM.index == wildcards.sample, 'bam'] + ".bai"
	output:
		org_header = temp(bam_dir + "{sample}.org.header"),
		new_header = temp(bam_dir + "{sample}.new.header"),
		bam = bam_dir + "{sample}.ready.bam",
		bai = bam_dir + "{sample}.ready.bam.bai"
	params:
		script = HOME_DIR + "workflow/scripts/reheader.py"
	threads: THREADS
	conda: NGS_env_1
	log: logs_dir + "link_bam/{sample}.log"
	priority: 0
	wildcard_constraints:
		sample = "|".join(bam_samples) if bam_samples else "no_wildcards"
	# run:
	# 	if not os.path.islink(output[0]):
	# 		os.symlink(input[0], output[0])
	# 	if not os.path.islink(output[1]):
	# 		os.symlink(input[1], output[1])
	shell:
		"(samtools view -H {input.bam} > {output.org_header}; "
		"python {params.script} {output.org_header} {output.new_header}; "
		"samtools reheader {output.new_header} {input.bam} > {output.bam}; "
		"samtools index -@ {threads} {output.bam}; "
		") > {log} 2>&1"



rule link_cram:
	input:
		cram = lambda wildcards: DF_BAM.loc[DF_BAM.index == wildcards.sample, 'bam'],
		crai = lambda wildcards: DF_BAM.loc[DF_BAM.index == wildcards.sample, 'bam'] + ".crai"
	output:
		org_header = temp(bam_dir + "{sample}.org.header"),
		new_header = temp(bam_dir + "{sample}.new.header"),
		cram = temp(bam_dir + "{sample}.ready.cram"),
		crai = temp(bam_dir + "{sample}.ready.cram.crai")
	params:
		script = HOME_DIR + "workflow/scripts/reheader.py"
	threads: THREADS
	conda: NGS_env_1
	log: logs_dir + "link_cram/{sample}.log"
	priority: 0
	wildcard_constraints:
		sample = "|".join(list(cram_samples)) if cram_samples else "no_wildcards"
	shell:
		# if not os.path.islink(output[0]):
		# 	os.symlink(input[0], output[0])
		# if not os.path.islink(output[1]):
		# 	os.symlink(input[1], output[1])
		"(samtools view -H {input.cram} > {output.org_header}; "
		"python {params.script} {output.org_header} {output.new_header}; "
		"samtools reheader {output.new_header} {input.cram} > {output.cram}; "
		"samtools index -@ {threads} {output.cram}; "
		") > {log} 2>&1"


rule cram2bam:
	input:
		cram = bam_dir + "{sample}.ready.cram",
		crai = bam_dir + "{sample}.ready.cram.crai"
	output:
		bam = bam_dir + "{sample}.ready.bam",
		bai = bam_dir + "{sample}.ready.bam.bai"
	params:
		ref = REFERENCE
	conda: NGS_env_1
	log: logs_dir + "cram2bam/{sample}.log"
	threads: THREADS
	priority: 0
	wildcard_constraints:
		sample = "|".join(cram_samples) if cram_samples else "no_wildcards"
	shell:
		"(samtools view -Sb -T {params.ref} -@ {threads} -o {output.bam} {input.cram}; "
		"samtools index -@ {threads} {output.bam} {output.bai}; "
		") > {log} 2>&1"
