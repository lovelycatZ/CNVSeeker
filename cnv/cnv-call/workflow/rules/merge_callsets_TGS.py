info_dir = merged_vcf_dir + "info/"
output_dir = merged_vcf_dir + "output/"


SAMPLES_PB = DF[DF['platform'].str.contains("pb")].index.tolist()
SAMPLES_ONT = DF[DF['platform'].str.contains("ont")].index.tolist()


combo_dict = {
	"pb": {
		"DEL":{
			"SS": "(DeBreak,Sniffles,cuteSV)(2)+(NanoVar,SVision_pro)(1)", # yes
			"S": "(NanoVar,SVision_pro)(1)", # yes
			"M": "(DeBreak,NanoVar,Sniffles)(2)+(SVision_pro)(1)", # yes
			"L": "(DeBreak,NanoVar,SVision_pro)(2)" # yes
		},
		"DUP":{}	
	},
	"ont": {
		"DEL":{
			"SS": "(DeBreak,NanoVar,cuteSV)(2)+(SVision_pro,Sniffles)(2)", # yes
			"S": "(DeBreak,SVision_pro,Sniffles)(2)+(NanoVar,cuteSV)(2)", # yes
			"M": "(DeBreak,NanoVar,SVision_pro)(2)", # yes
			"L": "(DeBreak,NanoVar)(2)+(SVision_pro,Sniffles,cuteSV)(3)" # yes
		},
		"DUP":{}
	}
}


def get_combo(wilcards):
	if wilcards.sample in SAMPLES_PB:
		platform = "pb"
	elif wilcards.sample in SAMPLES_ONT:
		platform = "ont"
	else:
		raise ValueError(f"Sample {wilcards.sample} platform not recognized!")
	cnv_type = wilcards.cnv_type
	bin_tag = wilcards.bin_tag
	return combo_dict[platform][f'{cnv_type}'][f'{bin_tag}']


rule standardize_vcf:
	input:
		vcf = tools_vcf_dir + "{tool}/output/{sample}_{tool}.vcf"
	output:
		vcf_SS = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.SS.vcf",
		vcf_S = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.S.vcf",
		vcf_M = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.M.vcf",
		vcf_L = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.L.vcf"
	params:
		scirpt = "/gpfs/hpc/home/lijc/xiangxud/project/new_result/hd_wgs/cnvseeker/standardize_vcf.py"
	shell:
		"(python {params.scirpt} {input.vcf} {output.vcf_SS} {wildcards.tool} {wildcards.sample} {wildcards.cnv_type} 50 500 True M; "
		"python {params.scirpt} {input.vcf} {output.vcf_S} {wildcards.tool} {wildcards.sample} {wildcards.cnv_type} 500 5000 True M; "
		"python {params.scirpt} {input.vcf} {output.vcf_M} {wildcards.tool} {wildcards.sample} {wildcards.cnv_type} 5000 100000 True M; "
		"python {params.scirpt} {input.vcf} {output.vcf_L} {wildcards.tool} {wildcards.sample} {wildcards.cnv_type} 100000 10000000 True M)"


rule merge_respective:
	input:
		vcfs = expand(info_dir + "{{sample}}/raw/{{sample}}.{tool}.{{cnv_type}}.{{bin_tag}}.vcf", tool = TOOLS)
	output:
		merged_vcf = info_dir + "{sample}/{sample}.{cnv_type}.{bin_tag}.merged.vcf",
		vcf_list = info_dir + "{sample}/{sample}.{cnv_type}.{bin_tag}.vcf_list",
		tmp_dir = directory(info_dir + "{sample}/merge_tmp.{sample}.{cnv_type}.{bin_tag}/")
	params:
		script = "/gpfs/hpc/home/lijc/xiangxud/project/new_result/hd_wgs/cnvseeker/suvivor_merge.py",
		# combo = lambda wildcards: combo_dict[NGS_DATA_TYPE][wildcards.cnv_type][wildcards.bin_tag],
		combo = lambda wildcards: get_combo(wildcards)
	conda: NGS_env_1
	threads: threads
	log: logs_dir + "merge_callsets/{sample}.{cnv_type}.{bin_tag}.log"
	shell:
		"(ls {input.vcfs} > {output.vcf_list}; "
		"python {params.script} {output.vcf_list} {output.merged_vcf} {output.tmp_dir} '{params.combo}'; "
		") > {log} 2>&1"


rule merge_final:
	input:
		vcfs = expand(info_dir + "{{sample}}/{{sample}}.{cnv_type}.{bin_tag}.merged.vcf", cnv_type = ["DEL"], bin_tag = ["SS", "S", "M", "L"])
	output:
		merged_vcf = info_dir + "{sample}/{sample}.merged.vcf",
		vcf_list = info_dir + "{sample}/{sample}.vcf_list",
	conda: NGS_env_1
	threads: threads
	log: logs_dir + "merge_callsets/{sample}.log"
	shell:
		"(ls {input.vcfs} > {output.vcf_list}; "
		"SURVIVOR merge {output.vcf_list} 0.5 1 1 1 0 50 {output.merged_vcf}; "
		") > {log} 2>&1"


rule add_id:
	input:
		filtered_vcf = info_dir + "{sample}/{sample}.merged.vcf"
	output:
		vcf = output_dir + "{sample}.vcf",
	run:
		add_var_id(input["filtered_vcf"], output["vcf"])
	

rule out_vcf:
	input:
		vcf = output_dir + "{sample}.vcf"
	output:
		vcf = final_dir + "{sample}.vcf"
	shell:	
		"cp {input.vcf} {output.vcf}"


rule out_bed:
	input:
		vcf = output_dir + "{sample}.vcf"
	output:
		bed = final_dir + "{sample}.bed"
	run:
		vcf2bed(input[0], output[0])
