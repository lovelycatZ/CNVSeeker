info_dir = merged_vcf_dir + "info/"
output_dir = merged_vcf_dir + "output/"



rule merge_1step:
	input:
		vcfs = expand(tools_vcf_dir + "{tool}/output/{{sample}}_{tool}.vcf", tool = TOOLS)
	output:
		#
		merged_vcf = info_dir + "{sample}/{sample}.merged.vcf",
		vcf_list = protected(info_dir + "{sample}/{sample}.vcf_list"),
		#
		filtered_vcf = info_dir + "{sample}/{sample}.merged.filtered.vcf",
		tmp1 = info_dir + "{sample}/{sample}.tmp1.vcf",
		# tmp2 = info_dir + "{sample}/{sample}.tmp2.vcf",	
		#	
		add_id_vcf = output_dir + "{sample}.vcf",
		#
		final_vcf = final_dir + "{sample}.vcf"
	params:
		script = HOME_DIR + "/workflow/scripts/suvivor_merge.1step.py",
		work_dir = info_dir + "{sample}/",
		data_type = DATA_TYPE.split("-")[-1],
		#
		script1 = HOME_DIR + "/workflow/scripts/blacklist_filter.py",
		script2 = HOME_DIR + "/workflow/scripts/gap_filter.py",
		# script3 = "/gpfs/hpc/home/lijc/xiangxud/project/new_result/hd_wgs/cnvseeker/segdup_filter.py",
		# script4 = "/gpfs/hpc/home/lijc/xiangxud/project/new_result/hd_wgs/cnvseeker/centromere_filter.py",
		s1 = config["base"]["black_list"]["value"],
		s2 = config["base"]["gap"]["value"],
		# s3 = "/gpfs/hpc/home/lijc/xiangxud/z_tmp_dataset/genomicSuperDups.b38.merge.bed",
		# s4 = "/gpfs/hpc/home/lijc/xiangxud/z_tmp_dataset/centromere.hg38.bed",
		script3 = HOME_DIR + "/workflow/scripts/add_id.py"
	conda: NGS_env_1
	threads: threads
	log: logs_dir + "merge_1step/{sample}.log"
	shell:
		#
		"(ls {input.vcfs} > {output.vcf_list}; "
		"python {params.script} {output.vcf_list} {output.merged_vcf} {params.work_dir} {params.data_type} {wildcards.sample}; "
		#
		"python {params.script1} -i {output.merged_vcf} -o {output.tmp1} -e {params.s1} -c 10 -p 0.7; "
		"python {params.script2} {output.tmp1} {params.s2} {output.filtered_vcf}; "
		#
		"python {params.script3} {output.filtered_vcf} {output.add_id_vcf}; "
		#
		"cp {output.add_id_vcf} {output.final_vcf}; "
		") > {log} 2>&1"
