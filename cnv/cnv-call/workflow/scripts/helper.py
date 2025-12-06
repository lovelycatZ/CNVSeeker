from ruamel import yaml
import pandas as pd
import pysam
import os
import re



def check_dir(dir):
	if not dir.endswith("/"):
		return dir + "/"
	else:
		return dir
	

	
def params_assembly(params_dict, concat_str):
	parmas_list = []
	for k, v in params_dict.items():
		if v is None:
			# print(f"value of {k} is set to None")
			continue
		
		if v is False:
			continue
		
		param = f"{k}"
		if not isinstance(v, bool):
			param += f"{concat_str}{v}"

		parmas_list.append(param)
	
	return " ".join(parmas_list)



def convert_config(dict_obj, top=None, key_chain_list=[]):
	"""自动替换字典中的模板字符串"""

	if top is None:
		top = dict_obj # top 保存原始的 config 文件

	for key, value in dict_obj.items():
		_key_chain_list = key_chain_list.copy() # 使用复制的列表进行后续处理，保证原列表由于引用传递不被修改

		if isinstance(value, dict):
			_key_chain_list.append(key)
			convert_config(value, top, _key_chain_list)
			
		if isinstance(value, str):
			# print('======= keys: ', key_chain_list, key)      # === DEBUG ===
			# print('old value: ', value)                       # === DEBUG ===
			tmpls = re.findall(r'\$\{[\.\w-]+?\}', value)
			for tmpl in tmpls:
				keys = tmpl[2:-1].split(".")
				val = top
				for k in keys:
					val = val[k]
				if type(val) != str and type(val) != yaml.scalarstring.DoubleQuotedScalarString and type(val) != yaml.scalarstring.SingleQuotedScalarString:
					print(f"(convert_config) Warning: The type of '{tmpl.strip('${} ')}' is not str, but will be forced to the str!")
					val = str(val)
				value = value.replace(tmpl, val)
			# print('new value: ', value)                       # === DEBUG ===

			update_dict = top # 用于 top 的更新
			for k in key_chain_list:
				update_dict = update_dict[k]
			# print('------- before：', '\n', top)              # === DEBUG ===
			update_dict[key] = value # 替换更新
			# print('------- after：', '\n', top)               # === DEBUG ===

		# TODO int, float, list, ...



def get_sm_tag(bam):
	with pysam.AlignmentFile(bam, "rb") as f:
		read_group = f.header.get("RG")
		if not read_group:
			raise Exception(f"Read group info not in file {bam}")
		read_group = read_group[0]
		if "SM" not in read_group:
			raise Exception(f"SM tag not in file {bam}")
		sm_tag = read_group["SM"]
		return sm_tag



def check_file(file):
	if not os.path.exists(file):
		raise Exception(f"File not exists: {file}, please check")
	


def check_platform(platform):
	if platform not in ["pb-hifi", "pb-clr", "ont"]:
		raise Exception(f"Platform not in correct format: {platform}, please check")
	


def check_sample_list_NGS(sample_list, tmp_dir):

	if not os.path.exists(tmp_dir):
		os.makedirs(tmp_dir)

	fq_list = tmp_dir + "fq.list"
	bam_list = tmp_dir + "bam.list"
	f_fq_list = open(fq_list, "w")
	f_bam_list = open(bam_list, "w")

	with open(sample_list, "r") as ff:
		for line in ff:
			# 跳过空行和以 # 开头的注释行
			if line.isspace() or line.startswith("#") :
				continue
			line_lst = line.rstrip().split("\t")
			sample_name = line_lst[0]
			if len(line_lst) == 4:
				f_bam_list.write(line)
				file = line_lst[3]
				if file.endswith(".bam"):
					file_index = file + ".bai"
				else:
					file_index = file + ".crai"	
				for f in [file, file_index]:
					check_file(f)
				sm_tag = get_sm_tag(file) 
				if sm_tag != sample_name:
					print(f"sample name for {file} was replaced with SM tag in BAM header as they were not consistent")
			else:
				f_fq_list.write(line)
				for f in line_lst[3:5]:
					check_file(f)
	f_bam_list.close()
	f_fq_list.close()



def parse_sample_list_NGS(sample_list, tmp_dir):

	check_sample_list_NGS(sample_list, tmp_dir)

	fq_list = tmp_dir + "fq.list"
	bam_list = tmp_dir + "bam.list"

	df_bam = pd.read_csv(bam_list, 
					sep = '\t', 
					comment = "#",
					names = ['sample', 'sex', 'tag', 'bam'], 
					index_col = 0)
	
	df_fq = pd.read_csv(fq_list, 
					sep = '\t', 
					comment = "#",
					names = ['sample', 'sex', 'tag', 'fq_r1', 'fq_r2'], 
					index_col = 0)

	df_sex = pd.concat([df_bam.iloc[:, :1], df_fq.iloc[:, :1]], axis = 0)
	
	return df_bam, df_fq, df_sex



def check_sample_list_TGS(sample_list):

	with open(sample_list, "r") as ff:
		for line in ff:
			# 跳过空行和以 # 开头的注释行
			if line.isspace() or line.startswith("#") :
				continue
			line_lst = line.rstrip().split("\t")
			sample_name, sex, platform, file = line_lst
			check_platform(platform)
			check_file(file)
			if file.endswith((".bam", ".cram")) and not file.endswith(".ccs.bam"):
				if file.endswith(".bam"):
					file_index = file + ".bai"
				else:
					file_index = file + ".crai"
				check_file(file_index)
				sm_tag = get_sm_tag(file)
				if sm_tag != sample_name:
					print(f"Sample name for {file} was replaced with SM tag in BAM header as they were not consistent")
			if file.endswith(".ccs.bam"):
				file_index = file + ".pbi"
				check_file(file_index)



def parse_sample_list_TGS(sample_list):
	
	check_sample_list_TGS(sample_list)

	df = pd.read_csv(sample_list, 
				  sep = '\t', 
				  names = ['sample', 'sex', 'platform', 'file'], 
				  comment = '#',
				  index_col = 0)
	return df



def make_dir_dict(tools, tools_vcf_dir):
	TOOLS_DIR_DICT = {}
	for tool in tools:
		TOOLS_DIR_DICT[tool] = {"output": tools_vcf_dir + tool  + "/output/", 
							"info": tools_vcf_dir + tool  + "/info/"}
	return TOOLS_DIR_DICT



def vcf2bed(vcf, bed):
	with open(bed, "w") as out:
		for line in open(vcf, "r"):
			if line.startswith("#"):
				continue
			line = line.rstrip().split("\t")
			sv_type = line[7].split("SVTYPE=")[1].split(";")[0]

			if sv_type not in ["DUP", "DEL"]:
				continue

			chr = line[0]
			if not chr.startswith("chr"):
				chr = "chr" + chr
				
			st = int(line[1]) - 1
			if line[7].startswith("END="):
				fl = line[7].split("END=")[1].split(";")[0]
			elif ";END=" in line[7]:
				fl = line[7].split(";END=")[1].split(";")[0]
			elif "SVLEN=" in line[7]:
				sv_len = int(line[7].split("SVLEN=")[1].split(";")[0])
				fl = st + sv_len
			else:
				continue

			print("\t".join([chr, str(st), str(fl), sv_type]), file = out)
