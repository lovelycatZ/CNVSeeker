## 质控
rule qc:
	input:
		read1 = lambda wildcards: get_fastq_sample(wildcards, "R1"),
		read2 = lambda wildcards: get_fastq_sample(wildcards, "R2")
	output:
		read1_qced = fastq_qced_dir + "{sample}_R1.fq.gz",
		read2_qced = fastq_qced_dir + "{sample}_R2.fq.gz",
		json_report = fastq_qced_dir + "{sample}.json",
		html_report = fastq_qced_dir + "{sample}.html"
	params:
		params_str = params_assembly(config["tools"]["fastp"]["params"], " ")
	conda: NGS_env_1
	benchmark: benchmark_dir + "qc/{sample}.log"
	log: logs_dir + "qc/{sample}.log"
	threads: THREADS
	shell:
		"(fastp -i {input.read1} -I {input.read2} "
		"-o {output.read1_qced} -O {output.read2_qced} "
		"-j {output.json_report} -h {output.html_report} "
		"{params.params_str}; "
		") > {log} 2>&1"
	

## 比对
rule align:
	input:
		read1_qced = fastq_qced_dir + "{sample}_R1.fq.gz",
		read2_qced = fastq_qced_dir + "{sample}_R2.fq.gz"
	output:
		sam = temp(bam_dir + "{sample}.sam")
	params:	
		ref = REFERENCE,
		rg = lambda wildcards: '"@RG\\tID:{sample}\\tPL:ILLUMINA\\tLB:{sample}\\tSM:{sample}"'.format(sample = wildcards.sample),
		params_str = params_assembly(config["tools"]["bwa"]["params"], " ")
	conda: NGS_env_1
	benchmark: benchmark_dir + "align/{sample}.log"
	log: logs_dir + "align/{sample}.log"
	threads: THREADS
	shell:
		"(bwa mem -R {params.rg} {params.params_str} {params.ref} {input.read1_qced} {input.read2_qced} > {output.sam}"
		") > {log} 2>&1"


# 转换为BAM格式
rule sam2bam:
	input:
		sam = bam_dir + "{sample}.sam"
	output:
		bam = temp(bam_dir + "{sample}.unsort.bam")
	conda: NGS_env_1
	benchmark: benchmark_dir + "sam2bam/{sample}.log"
	log: logs_dir + "sam2bam/{sample}.log"
	threads: THREADS
	shell:
		"(samtools view -bS -@ {threads} {input.sam} > {output.bam};"
		") > {log} 2>&1"


# 排序
rule sort:
	input:
		bam = rules.sam2bam.output.bam
	output:
		bam = temp(bam_dir + "{sample}.sort.bam"),
		samtools_sort_tmp_dir = temp(directory(bam_dir + "samtools_sort_tmp_{sample}"))
	conda: NGS_env_1
	benchmark: benchmark_dir + "sort/{sample}.log"
	log: logs_dir + "sort/{sample}.log"
	threads: THREADS
	shell:
		"(mkdir -p {output.samtools_sort_tmp_dir} && cd {output.samtools_sort_tmp_dir}; "
		"samtools sort -@ {threads} {input.bam} > {output.bam}; "
		") > {log} 2>&1"


## 标记(去除)重复序列
rule markdup:
	input:
		bam = rules.sort.output.bam
	output:
		bam = temp(bam_dir + "{sample}.sort.mkdup.bam")
	params:
		params_str = params_assembly(config["tools"]["sambamba"]["params"], " ")
	conda: NGS_env_1
	benchmark: benchmark_dir + "markdup/{sample}.log"
	log: logs_dir + "markdup/{sample}.log"
	threads: THREADS
	shell:
		"(sambamba markdup {params.params_str} {input.bam} {output.bam}; "
		") > {log} 2>&1"


# BaseRecalibrator
rule BQSR_create:
	input:
		bam = rules.markdup.output.bam
	output:
		table = bam_dir + "{sample}.BQSR.table"
	params:
		ref = REFERENCE,
		known_site_1 = config["tools"]["gatk4"]["params"][f"--known_sites_1_{BUILD}"],
		known_site_2 = config["tools"]["gatk4"]["params"][f"--known_sites_2_{BUILD}"]
	conda: gatk4
	benchmark: benchmark_dir + "BQSR_create/{sample}.log"
	log: logs_dir + "BQSR_create/{sample}.log"
	threads: THREADS
	shell:
		"(gatk BaseRecalibrator -R {params.ref} -I {input.bam} -O {output.table} "
		"--known-sites {params.known_site_1} "
		"--known-sites {params.known_site_2}; "
		") > {log} 2>&1"


# ApplyBQSR
rule BQSR_apply:
	input:
		bam = rules.markdup.output.bam,
		table = rules.BQSR_create.output.table
	output:
		bam = bam_dir + "{sample}.ready.bam",
		bai = bam_dir + "{sample}.ready.bam.bai"
	params:
		ref = REFERENCE,
		bai = bam_dir + "{sample}.ready.bai"
	conda: gatk4
	benchmark: benchmark_dir + "BQSR_apply/{sample}.log"
	log: logs_dir + "BQSR_apply/{sample}.log"
	threads: THREADS
	wildcard_constraints:
		sample = "|".join(fq_samples + fq_gz_samples) if fq_samples + fq_gz_samples else "no_wildcards"
	shell:
		"(gatk ApplyBQSR -R {params.ref} -I {input.bam} -O {output.bam} --bqsr-recal-file {input.table}; "
		"mv {params.bai} {output.bai}; "
		") > {log} 2>&1"


rule bam2cram:
	input:
		bam = bam_dir + "{sample}.ready.bam",
		bai = bam_dir + "{sample}.ready.bam.bai"
	output:
		cram = bam_dir + "{sample}.ready.cram",
		crai = bam_dir + "{sample}.ready.cram.crai"
	params:
		ref = REFERENCE
	conda: NGS_env_1
	benchmark: benchmark_dir + "cram2bam/{sample}.log"
	log: logs_dir + "cram2bam/{sample}.log"
	threads: THREADS
	wildcard_constraints:
		sample = "|".join(fq_samples + fq_gz_samples) if fq_samples + fq_gz_samples else "no_wildcards"
	shell:
		"(samtools view -SC -T {params.ref} -@ {threads} {input[0]} > {output.cram}; "
		"samtools index -@ {threads} {output.cram} {output.crai}; "
		") > {log} 2>&1"


rule out_cram:
	input:
		cram = bam_dir + "{sample}.ready.cram",
		crai = bam_dir + "{sample}.ready.cram.crai"
	output:
		cram = final_dir + "{sample}.cram",
		crai = final_dir + "{sample}.cram.crai"
	benchmark: benchmark_dir + "out_cram/{sample}.log"
	wildcard_constraints:
		sample = "|".join(fq_samples + fq_gz_samples) if fq_samples + fq_gz_samples else "no_wildcards"
	shell:	
		# "(cp {input.cram} {output.cram}; "
		# "cp {input.crai} {output.crai})"
		"(ln -s {input.cram} {output.cram}; "
		"ln -s {input.crai} {output.crai})"


rule out_bam:
	input:
		bam = bam_dir + "{sample}.ready.bam",
		bai = bam_dir + "{sample}.ready.bam.bai"
	output:
		bam = final_dir + "{sample}.bam",
		bai = final_dir + "{sample}.bam.bai"
	benchmark: benchmark_dir + "out_bam/{sample}.log"
	shell:	
		# "(cp {input.bam} {output.bam}; "
		# "cp {input.bai} {output.bai})"
		"(ln -s {input.bam} {output.bam}; "
		"ln -s {input.bai} {output.bai})"
