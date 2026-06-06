info_dir = merged_vcf_dir + "info/"
output_dir = merged_vcf_dir + "output/"


COMBO_DICT = {
	"HD_WGS":{
		"DEL":{
			"SS": "(Wham,delly)(2)+(lumpy)(1)", # yes
			"S": "(delly,lumpy)(2)", # yes
			"M": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)", # yes
            "L": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)"
		},
		"DUP":{
			"SS": "(Wham)(1)+(delly,lumpy)(2)", # yes
			"S": "(Wham,lumpy)(2)+(delly)(1)", # yes
			"M": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)", # yes
            "L": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)"
		}
	},
	"LD_WGS":{
		"DEL":{
			"SS": "(delly,gridss,lumpy)(1)", # yes
			"S": "(delly,gridss,lumpy)(1)", # yes
			"M": "(CNVpytor,cn.mops,delly,lumpy)(2)+(gridss)(1)", # yes
			"L": "(CNVpytor,delly)(2)+(ControlFREEC,cn.mops)(2)" # yes
		},
		"DUP":{
			"SS": "(Manta,Wham,gridss)(1)", # yes
			"S": "(Manta,gridss)(1)", # yes
			"M": "(CNVpytor,Manta,Wham)(1)+(ControlFREEC,gridss)(1)", # yes
			"L": "(CNVpytor)(1)+(ControlFREEC,Manta,Wham,cn.mops)(3)", # yes
		}
	},
	"WES":{
		"DEL":{
			"S": "(ECOLE,ExomeDepth,gatk4,xhmm)(2)+(cn.mops)(1)", # yes
			"M": "(ECOLE,ExomeDepth,cn.mops,gatk4,xhmm)(1)", # yes
			"L": "(ECOLE,gatk4)(1)+(ExomeDepth,xhmm)(1)" # yes
		},
		"DUP":{
			"S": "(ECOLE)(1)+(ExomeDepth)(1)", # yes
			"M": "(ECOLE)(1)+(ExomeDepth,gatk4,xhmm)(2)", # yes
			"L": "(ECOLE,ExomeDepth,cn.mops,xhmm)(1)", # yes
		}
	}
}


if DATA_TYPE == "NGS-WES":
	bin_tags = ["S", "M", "L"]
else:
	bin_tags = ['SS', 'S', 'M', 'L']

cnv_types = ['DEL', 'DUP']



if DATA_TYPE == "NGS-WES":
	rule standardize_vcf:
		input:
			vcf = tools_vcf_dir + "{tool}/output/{sample}_{tool}.vcf"
		output:
			vcf_S = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.S.vcf",
			vcf_M = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.M.vcf",
			vcf_L = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.L.vcf"
		params:
			scirpt = HOME_DIR + "/workflow/scripts/standardize_vcf.py"
		shell:
			"(python {params.scirpt} {input.vcf} {output.vcf_S} {wildcards.tool} {wildcards.sample} {wildcards.cnv_type} 50 10000 True M; "
			"python {params.scirpt} {input.vcf} {output.vcf_M} {wildcards.tool} {wildcards.sample} {wildcards.cnv_type} 10000 100000 True M; "
			"python {params.scirpt} {input.vcf} {output.vcf_L} {wildcards.tool} {wildcards.sample} {wildcards.cnv_type} 100000 10000000 True M)"
else:
	rule standardize_vcf:
		input:
			vcf = tools_vcf_dir + "{tool}/output/{sample}_{tool}.vcf"
		output:
			vcf_SS = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.SS.vcf",
			vcf_S = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.S.vcf",
			vcf_M = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.M.vcf",
			vcf_L = info_dir + "{sample}/raw/{sample}.{tool}.{cnv_type}.L.vcf"
		params:
			scirpt = HOME_DIR + "/workflow/scripts/standardize_vcf.py"
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
		script = HOME_DIR + "/workflow/scripts/suvivor_merge.py",
		combo = lambda wildcards: COMBO_DICT[NGS_DATA_TYPE][wildcards.cnv_type][wildcards.bin_tag]
	conda: NGS_env_1
	threads: THREADS
	log: logs_dir + "merge_callsets/{sample}.{cnv_type}.{bin_tag}.log"
	shell:
		"(ls {input.vcfs} > {output.vcf_list}; "
		"python {params.script} {output.vcf_list} {output.merged_vcf} {output.tmp_dir} '{params.combo}'; "
		") > {log} 2>&1"


rule merge_final:
	input:
		vcfs = expand(info_dir + "{{sample}}/{{sample}}.{cnv_type}.{bin_tag}.merged.vcf", cnv_type = cnv_types, bin_tag = bin_tags)
	output:
		merged_vcf = info_dir + "{sample}/{sample}.merged.vcf",
		vcf_list = info_dir + "{sample}/{sample}.vcf_list",
	conda: NGS_env_1
	threads: THREADS
	log: logs_dir + "merge_callsets/{sample}.log"
	shell:
		"(ls {input.vcfs} > {output.vcf_list}; "
		"SURVIVOR merge {output.vcf_list} 0.5 1 1 1 0 50 {output.merged_vcf}; "
		") > {log} 2>&1"


rule post_filter:
	input:
		merged_vcf = info_dir + "{sample}/{sample}.merged.vcf"
	output:
		filtered_vcf = info_dir + "{sample}/{sample}.merged.filtered.vcf",
		tmp1 = info_dir + "{sample}/{sample}.tmp1.vcf",
		# tmp2 = info_dir + "{sample}/{sample}.tmp2.vcf",
		# tmp3 = info_dir + "{sample}/{sample}.tmp3.vcf",
	params:
		script1 = HOME_DIR + "/workflow/scripts/blacklist_filter.py",
		script2 = HOME_DIR + "/workflow/scripts/gap_filter.py",
		# script3 = HOME_DIR + "/workflow/scripts/segdup_filter.py",
		# script4 = HOME_DIR + "/workflow/scripts/centromere_filter.py",
		s1 = config["base"]["black_list"]["value"],
		s2 = config["base"]["gap"]["value"],
		# s3 = "/gpfs/hpc/home/lijc/xiangxud/z_tmp_dataset/genomicSuperDups.b38.merge.bed",
		# s4 = "/gpfs/hpc/home/lijc/xiangxud/z_tmp_dataset/centromere.hg38.bed",
	conda: NGS_env_1
	threads: THREADS
	log: logs_dir + "post_filter/{sample}.log"
	shell:
		"(python {params.script1} -i {input.merged_vcf} -o {output.tmp1} -e {params.s1} -c 10 -p 0.7; "
		"python {params.script2} {output.tmp1} {params.s2} {output.filtered_vcf}; "
		# "python {params.script3} -i {output.tmp2} -o {output.filtered_vcf} -e {params.s3} -c 10 -p 0.7; "
		# "python {params.script4} {params.s4} {output.tmp3} {output.filtered_vcf}; "
		") > {log} 2>&1"


rule add_id:
	input:
		filtered_vcf = info_dir + "{sample}/{sample}.merged.filtered.vcf"
	output:
		vcf = output_dir + "{sample}.vcf",
	params:
		script = HOME_DIR + "/workflow/scripts/add_id.py"
	shell:
		"python {params.script3} {output.filtered_vcf} {output.add_id_vcf}; "
	

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
