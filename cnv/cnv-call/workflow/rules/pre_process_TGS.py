# rule qc:
# 	input:
# 		fastq = lambda wildcards: get_fastq_sample(wildcards)
# 	output:
# 		fastq_qced = fastq_qced_dir + "{sample}.{platform}.fq.gz"
# 	params:
# 		params_str = params_assembly(config["tools"]["filtlong"]["params"], " ")
# 	conda: TGS_env_1
# 	log: logs_dir + "filtlong_qc/{sample}.{platform}.log"
# 	threads: THREADS
# 	wildcard_constraints:
# 		sample = "|".join(fq_samples + fq_gz_samples) if fq_samples + fq_gz_samples else "no_wildcards"
# 	shell:
# 		"(filtlong {params.params_str} {input.fastq} | gzip > {output.fastq_qced}; "
# 		") > {log} 2>&1"


rule minimap2_align:
	input:
		fastq = lambda wildcards: get_fastq_sample(wildcards)
	output:
		sam = temp(bam_dir + "{sample}.{platform}.minimap2.sam")
	params:
		ref = REFERENCE,
		pt_param = lambda wildcards: get_aligner_params(wildcards, "minimap2"),
		rg = lambda wildcards: '"@RG\\tID:{sample}\\tSM:{sample}"'.format(sample = wildcards.sample)
	conda: TGS_env_1
	benchmark: benchmark_dir + "minimap2_align/{sample}.{platform}.log"
	log: logs_dir + "minimap2_align/{sample}.{platform}.log"
	threads: THREADS
	# wildcard_constraints:
	# 	sample = "|".join(fq_samples + fq_gz_samples) if fq_samples + fq_gz_samples else "no_wildcards"
	shell:
		"(minimap2 {params.pt_param} --MD -Y -a -t {threads} -R {params.rg} {params.ref} {input.fastq} > {output.sam}; "
		") > {log} 2>&1"


rule minimap2_align_asm:
	input:
		assmb_fasta = lambda wildcards: get_fasta_sample(wildcards)
	output:
		sam = temp(bam_dir + "{sample}.{platform}.minimap2-asm.sam")
	params:
		ref = REFERENCE
	conda: TGS_env_1
	benchmark: benchmark_dir + "minimap2_align/{sample}.{platform}.log"
	log: logs_dir + "minimap2_align/{sample}.{platform}.log"
	threads: THREADS
	# wildcard_constraints:
	# 	sample = "|".join(fa_samples + fa_gz_samples) if fa_samples + fa_gz_samples else "no_wildcards"
	shell:
		"(minimap2 -a --cs -r2k --eqx -x asm5 -t {threads} {params.ref} {input.assmb_fasta} > {output.sam}; "
		") > {log} 2>&1"


# rule ngmlr_align:
# 	input:
# 		fastq = fastq_qced_dir + "{sample}.{platform}.fq.gz"
# 	output:
# 		sam = temp(bam_dir + "{sample}.{platform}.ngmlr.sam")
# 	params:
# 		ref = REFERENCE,
# 		pt_param = lambda wildcards: get_aligner_params(wildcards, "ngmlr"),
# 		rg = lambda wildcards: "--rg-id {sample} --rg-sm {sample}".format(sample = wildcards.sample)
# 	conda: TGS_env_1
# 	log: logs_dir + "ngmlr_align/{sample}.{platform}.log"
# 	threads: THREADS
# 	shell:
# 		"(ngmlr {params.rg} {params.pt_param} -t {threads} -r {params.ref} -q {input.fastq} -o {output.sam}; "
# 		") > {log} 2>&1"


rule pbmm2_align:
	input:
		# fastq = fastq_qced_dir + "{sample}.{platform}.fq.gz"
		fastq = lambda wildcards: get_fastq_sample(wildcards)
	output:
		bam = bam_dir + "{sample}.{platform}.pbmm2.sort.bam",
		bai = bam_dir + "{sample}.{platform}.pbmm2.sort.bam.bai"
	params:
		ref = REFERENCE,
		pt_param = lambda wildcards: get_aligner_params(wildcards, "pbmm2"),
		rg = lambda wildcards: '"@RG\\tID:{sample}\\tSM:{sample}"'.format(sample = wildcards.sample)
	conda: TGS_env_1
	benchmark: benchmark_dir + "pbmm2_align/{sample}.{platform}.log"
	log: logs_dir + "pbmm2_align/{sample}.{platform}.log"
	threads: THREADS
	wildcard_constraints:
		platform = "|".join(["pb-hifi", "pb-clr"])
	shell:
		"(pbmm2 align {params.pt_param} --rg {params.rg} --sort -j {threads} -J {threads} {params.ref} {input.fastq} {output.bam}; "
		") > {log} 2>&1"


# rule meryl_kmers:
# 	input:
# 		ref = REFERENCE
# 	output:
# 		merylDB = directory(bam_dir + "merylDB"),
# 		kmers = bam_dir + "kmers.txt"
# 	conda: TGS_env_1
# 	log: logs_dir + "winnowmap_align/meryl_kmers.log"
# 	threads: THREADS
# 	shell:
# 		"(meryl count k=15 output {output.merylDB} {input.ref}; "
# 		"meryl print greater-than distinct=0.9998 {output.merylDB} > {output.kmers}; "
# 		") > {log} 2>&1"


# rule winnowmap_align:
# 	input:
# 		fastq = fastq_qced_dir + "{sample}.{platform}.fq.gz",
# 		kmers = bam_dir + "kmers.txt"
# 	output:
# 		sam = temp(bam_dir + "{sample}.{platform}.winnowmap.sam")
# 	params:
# 		ref = REFERENCE,
# 		pt_param = lambda wildcards: get_aligner_params(wildcards, "winnowmap"),
# 		rg = lambda wildcards: '"@RG\\tID:{sample}\\tSM:{sample}"'.format(sample = wildcards.sample)
# 	conda: TGS_env_1
# 	log: logs_dir + "winnowmap_align/{sample}.{platform}.log"
# 	threads: THREADS
# 	shell:
# 		"(winnowmap {params.pt_param} -W {input.kmers} -a -t {threads} {params.ref} {input.fastq} > {output.sam}; "
# 		") > {log} 2>&1"


rule sam2bam:
	input:
		sam = bam_dir + "{sample}.{platform}.{aligner}.sam"
	output:
		bam = temp(bam_dir + "{sample}.{platform}.{aligner}.bam")
	conda: TGS_env_1
	benchmark: benchmark_dir + "sam2bam/{sample}.{platform}.{aligner}.log"
	log: logs_dir + "sam2bam/{sample}.{platform}.{aligner}.log"
	threads: THREADS
	shell:
		"(samtools view -bS -@ {threads} {input.sam} > {output.bam}; "
		") > {log} 2>&1"


rule sort_bam:
	input:
		bam = bam_dir + "{sample}.{platform}.{aligner}.bam"
	output:	
		bam = bam_dir + "{sample}.{platform}.{aligner}.sort.bam",
		bai = bam_dir + "{sample}.{platform}.{aligner}.sort.bam.bai"
	conda: TGS_env_1
	benchmark: benchmark_dir + "sort_bam/{sample}.{platform}.{aligner}.log"
	log: logs_dir + "sort_bam/{sample}.{platform}.{aligner}.log"
	threads: THREADS
	shell:
		"(samtools sort -@ {threads} {input.bam} -o {output.bam}; "
		"samtools index -@ {threads} {output.bam}; "
		") > {log} 2>&1"


rule archive_ready_bam:
	input:
		bam = bam_dir + "{sample}.{platform}.{aligner}.sort.bam",
		bai = bam_dir + "{sample}.{platform}.{aligner}.sort.bam.bai"
	output:
		bam = bam_dir + "{sample}.{platform}.{aligner}.ready.bam",
		bai = bam_dir + "{sample}.{platform}.{aligner}.ready.bam.bai"
	threads: THREADS
	benchmark: benchmark_dir + "archive_ready_bam/{sample}.{platform}.{aligner}.log"
	log: logs_dir + "archive_ready_bam/{sample}.{platform}.{aligner}.log"
	wildcard_constraints:
		sample = "|".join(fa_samples + fa_gz_samples + fq_gz_samples + fq_samples + ccs_bam_samples) if fa_samples + fa_gz_samples + fq_gz_samples + fq_samples + ccs_bam_samples else "no_wildcards"
	run:
		if not os.path.islink(output[0]):
			os.symlink(input[0], output[0])
		if not os.path.islink(output[1]):
			os.symlink(input[1], output[1])
		

rule bam2cram:
	input:
		bam = bam_dir + "{sample}.{platform}.{aligner}.sort.bam",
		bai = bam_dir + "{sample}.{platform}.{aligner}.sort.bam.bai"
	output:
		cram = bam_dir + "{sample}.{platform}.{aligner}.ready.cram",
		crai = bam_dir + "{sample}.{platform}.{aligner}.ready.cram.crai"
	params:
		ref = REFERENCE
	conda: TGS_env_1
	benchmark: benchmark_dir + "cram2bam/{sample}.{platform}.{aligner}.log"
	log: logs_dir + "cram2bam/{sample}.{platform}.{aligner}.log"
	threads: THREADS
	shell:
		"(samtools view -SC -T {params.ref} -@ {threads} {input[0]} > {output.cram}; "
		"samtools index -@ {threads} {output.cram} {output.crai}; "
		") > {log} 2>&1"


rule out_cram_1:
	input:
		cram = bam_dir + "{sample}.{platform}.ready.cram",
		crai = bam_dir + "{sample}.{platform}.ready.cram.crai"
	output:
		cram = final_dir + "{sample}.{platform}.cram",
		crai = final_dir + "{sample}.{platform}.cram.crai"
	benchmark: benchmark_dir + "out_cram_1/{sample}.{platform}.log"
	shell:	
		"(cp {input.cram} {output.cram}; "
		"cp {input.crai} {output.crai})"


rule out_cram_2:
	input:
		cram = bam_dir + "{sample}.{platform}.{aligner}.ready.cram",
		crai = bam_dir + "{sample}.{platform}.{aligner}.ready.cram.crai"
	output:
		cram = final_dir + "{sample}.{platform}.{aligner}.cram",
		crai = final_dir + "{sample}.{platform}.{aligner}.cram.crai"
	benchmark: benchmark_dir + "out_cram_2/{sample}.{platform}.{aligner}.log"
	shell:	
		"(cp {input.cram} {output.cram}; "
		"cp {input.crai} {output.crai})"
		

rule out_bam_1:
	input:
		bam = bam_dir + "{sample}.{platform}.ready.bam",
		bai = bam_dir + "{sample}.{platform}.ready.bam.bai"
	output:
		bam = final_dir + "{sample}.{platform}.bam",
		bai = final_dir + "{sample}.{platform}.bam.bai"
	benchmark: benchmark_dir + "out_bam_1/{sample}.{platform}.log"
	shell:	
		"(cp {input.bam} {output.bam}; "
		"cp {input.bai} {output.bai})"


rule out_bam_2:
	input:
		bam = bam_dir + "{sample}.{platform}.{aligner}.ready.bam",
		bai = bam_dir + "{sample}.{platform}.{aligner}.ready.bam.bai"
	output:
		bam = final_dir + "{sample}.{platform}.{aligner}.bam",
		bai = final_dir + "{sample}.{platform}.{aligner}.bam.bai"
	benchmark: benchmark_dir + "out_bam_2/{sample}.{platform}.{aligner}.log"
	shell:	
		"(cp {input.bam} {output.bam}; "
		"cp {input.bai} {output.bai})"


# rule subreads2ccs:
# 	input:
# 		subreads_bam = subreads_dir + "{sample}.subreads.bam"
# 	output:
# 		ccs_bam = ccs_dir + "{sample}.ccs.bam"
# 	conda: TGS_env_1
# 	log: bam_dir + "subreads2ccs/{sample}.log"
# 	threads: THREADS
# 	shell:
# 		"(ccs --hifi-kinetics -j {threads} {input.subreads_bam} {output.ccs_bam}; "
# 		") > {log} 2>&1"
