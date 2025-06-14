info_dir = merged_vcf_dir + "info/"
output_dir = merged_vcf_dir + "output/"


def get_vcfs(wildcards):
	vcf_list = []
	sample = wildcards.sample
	for tool in TOOLS:
		platfom = str(DF.loc[sample, 'platform'])
		if tool == "pbsv" and platfom == "ont":
			continue
		vcf = info_dir + f"{sample}/raw/{sample}_{tool}.vcf"
		vcf_list.append(vcf)
	return vcf_list


for sample in SAMPLES:
	platfom = str(DF.loc[sample, 'platform'])
	DF_tool = DF.copy(deep=True)
	rule:
		input:
			vcf = tools_vcf_dir + "{{tool}}/output/{sample}_{{tool}}.vcf".format(sample = sample)
		output:
			vcf = info_dir + "{sample}/raw/{sample}_{{tool}}.vcf".format(sample = sample)
		params:
			chrom_list = CHROM_LIST,
			sample_name = sample,
			tool = "{tool}"
		benchmark: benchmark_dir + "pre_merge/{sample}_{{tool}}.log".format(sample = sample)
		run:
			in_vcf = input["vcf"]
			out_vcf = output["vcf"]
			chrom_list = params["chrom_list"]
			sample_name = params["sample_name"]
			tool = params["tool"]
			sex = DF_tool.loc[sample_name, "sex"]
			make_pre_merge_vcf(in_vcf, out_vcf, sex, tool, chrom_list)


rule TGS_merge:
	input:
		# vcfs = expand(info_dir + "{{sample}}/raw/{{sample}}_{tool}.vcf", tool = TOOLS_1),
		vcfs = lambda wildcards: get_vcfs(wildcards)
	output:
		merged_vcf = info_dir + "{sample}/{sample}.merged.vcf",
	params:
		nanovar_vcf = info_dir + "{sample}/raw/{sample}_NanoVar.vcf",
		vcf_list = info_dir + "{sample}/{sample}.vcf_list",
		tmp_merged_vcf = info_dir + "{sample}/{sample}.merged.tmp.vcf",
		tmp_vcf_list = info_dir + "{sample}/{sample}.vcf_list_tmp",
		work_dir = info_dir + "{sample}/"
	conda: NGS_env_2
	benchmark: benchmark_dir + "TGS_merge/{sample}.log"
	threads: THREADS
	log: logs_dir + "TGS_merge/{sample}.log"
	shell:
		# "(ls {input.vcfs} > {params.tmp_vcf_list}; "
		# "SURVIVOR merge {params.tmp_vcf_list} 0.5 3 1 1 0 50 {params.tmp_merged_vcf}; "
		# "ls {params.nanovar_vcf} {params.tmp_merged_vcf} > {params.vcf_list}; "
		# "SURVIVOR merge {params.vcf_list} 0.1 1 1 1 0 50 {output.merged_vcf}; "
		# ") > {log} 2>&1"
		"(ls {input.vcfs} > {params.tmp_vcf_list}; "
		"jasmine file_list={params.tmp_vcf_list} out_file={params.tmp_merged_vcf} out_dir={params.work_dir} threads={threads} "
		"min_support=3 max_dist_linear=1 min_dist=50 min_overlap=0.5 --ignore_strand; "
		"ls {params.nanovar_vcf} {params.tmp_merged_vcf} > {params.vcf_list}; "
		"jasmine file_list={params.vcf_list} out_file={output.merged_vcf} out_dir={params.work_dir} threads={threads} "
		"min_support=1 max_dist_linear=1 min_dist=50 min_overlap=0.9 --ignore_strand; "
		") > {log} 2>&1"


rule post_filter:
	input:
		merged_vcf = info_dir + "{sample}/{sample}.merged.vcf"
	output:
		filtered_vcf = info_dir + "{sample}/{sample}.merged.filtered.vcf"
	params:
		script = HOME_DIR + "workflow/scripts/post_filter.py",
		black_list = config["base"]["black_list"]["value"],
		gap = config["base"]["gap"]["value"],
		work_dir = info_dir + "{sample}/filter_tmp/",
	conda: NGS_env_1
	benchmark: benchmark_dir + "post_filter/{sample}.log"
	threads: THREADS
	log: logs_dir + "post_filter/{sample}.log"
	shell:
		"(mkdir -p {params.work_dir}; cd {params.work_dir}; "
		"python {params.script} -e {params.black_list} -i {input.merged_vcf} -o tmp1.vcf -p 0.8; "
		"python {params.script} -e {params.gap} -i tmp1.vcf -o {output.filtered_vcf} -p 0.8; "
		") > {log} 2>&1"


rule add_id:
	input:
		vcf = info_dir + "{sample}/{sample}.merged.filtered.vcf"
	output:
		vcf = output_dir + "{sample}.vcf"
	benchmark: benchmark_dir + "add_id/{sample}.log"
	run:
		add_var_id(input["vcf"], output["vcf"])


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
