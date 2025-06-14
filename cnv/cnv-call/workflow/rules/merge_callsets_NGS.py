import os


info_dir = merged_vcf_dir + "info/"
output_dir = merged_vcf_dir + "output/"



# rule pre_merge_vcf:
# 	input:
# 		vcfs = expand(tools_vcf_dir + "{tool}/output/{{sample}}_{tool}.vcf", tool = TOOLS)
# 	output:
# 		vcfs = expand(info_dir + "{{sample}}/raw/{{sample}}_{tool}.vcf", tool = TOOLS)
# 	params:
# 		sample_name = "{sample}"
# 	run:
# 		input_list = input["vcfs"]
# 		output_list = output["vcfs"]
# 		sample_name = params["sample_name"]
# 		for in_vcf, out_vcf in zip(input_list, output_list):
# 			for tool in TOOLS:
# 				if tool in in_vcf:
# 					break
# 			make_pre_merge_vcf(in_vcf, out_vcf, sample_name, tool)



for tool in TOOLS:
	DF_tool = DF_SEX.copy(deep=True)
	rule:
		input:
			vcf = tools_vcf_dir + "{tool}/output/{{sample}}_{tool}.vcf".format(tool = tool)
		output:
			vcf = info_dir + "{{sample}}/raw/{{sample}}_{tool}.vcf".format(tool = tool)
		params:
			chrom_list = CHROM_LIST,
			sample_name = "{sample}",
			tool = tool
		benchmark: benchmark_dir + "pre_merge/{{sample}}_{tool}.log".format(tool = tool)
		run:
			in_vcf = input["vcf"]
			out_vcf = output["vcf"]
			chrom_list = params["chrom_list"]
			sample_name = params["sample_name"]
			tool = params["tool"]
			sex = DF_tool.loc[sample_name, "sex"]
			make_pre_merge_vcf(in_vcf, out_vcf, sex, tool, chrom_list)


rule NGS_merge:
	input:
		vcfs = expand(info_dir + "{{sample}}/raw/{{sample}}_{tool}.vcf", tool = TOOLS)
	output:
		vcf = info_dir + "{sample}/{sample}.merged.vcf"
	params:
		vcf_list = info_dir + "{sample}/vcf_list",
		work_dir = info_dir + "{sample}/"
	conda: NGS_env_2
	benchmark: benchmark_dir + "NGS_merge/{sample}.log"
	threads: THREADS
	log: logs_dir + "NGS_merge/{sample}.log"
	shell:
		"(ls {input.vcfs} > {params.vcf_list}; "
		"jasmine file_list={params.vcf_list} out_file={output.vcf} out_dir={params.work_dir} threads={threads} "
		"min_support=1 max_dist_linear=1 min_dist=50 min_overlap=0.5 --ignore_strand; "
		") > {log} 2>&1"
	

rule post_filter:
	input:
		merged_vcf = info_dir + "{sample}/{sample}.merged.vcf"
	output:
		filtered_vcf = info_dir + "{sample}/{sample}.merged.filtered.vcf"
	params:
		script1 = HOME_DIR + "workflow/scripts/merge_filter_wes.py" if NGS_DATA_TYPE == "WES" else HOME_DIR + "workflow/scripts/merge_filter_wgs.py",
		script2 = HOME_DIR + "workflow/scripts/post_filter.py",
		black_list = config["base"]["black_list"]["value"],
		gap = config["base"]["gap"]["value"],
		work_dir = info_dir + "{sample}/filter_tmp/",
		bed = EXOME_INTERVALS,
		ngs_data_type = NGS_DATA_TYPE,
	conda: NGS_env_1
	benchmark: benchmark_dir + "post_filter/{sample}.log"
	threads: THREADS
	log: logs_dir + "post_filter/{sample}.log"
	shell:
		"(mkdir -p {params.work_dir}; cd {params.work_dir}; "
		"python {params.script1} {input.merged_vcf} tmp1.vcf filtered.vcf {params.ngs_data_type} {params.bed}; "
		"python {params.script2} -e {params.black_list} -i tmp1.vcf -o tmp2.vcf -p 0.8; "
		"python {params.script2} -e {params.gap} -i tmp2.vcf -o {output.filtered_vcf} -p 0.8; "
		") > {log} 2>&1"


rule add_id:
	input:
		filtered_vcf = info_dir + "{sample}/{sample}.merged.filtered.vcf"
	output:
		vcf = output_dir + "{sample}.vcf",
	benchmark: benchmark_dir + "add_id/{sample}.log"
	run:
		add_var_id(input["filtered_vcf"], output["vcf"])
	

rule out_vcf:
	input:
		vcf = output_dir + "{sample}.vcf"
	output:
		vcf = final_dir + "{sample}.vcf"
	benchmark: benchmark_dir + "out_vcf/{sample}.log"
	shell:	
		"cp {input.vcf} {output.vcf}"


rule out_bed:
	input:
		vcf = output_dir + "{sample}.vcf"
	output:
		bed = final_dir + "{sample}.bed"
	benchmark: benchmark_dir + "out_bed/{sample}.log"
	run:
		vcf2bed(input[0], output[0])
