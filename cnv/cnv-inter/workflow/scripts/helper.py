def check_dir(dir):
	if not dir.endswith("/"):
		return dir + "/"
	else:
		return dir


def filt_wildcards(samples, suffixs):
	new_samples = []
	new_suffixs = []
	for sample, suffix in zip(samples, suffixs):
		if suffix in ["bed", "vcf"]:
			new_samples.append(sample)
			new_suffixs.append(suffix)
	return new_samples, new_suffixs


def vcf2bed(vcf, bed):

	out = open(bed, "w")
	print("\t".join(["#chr", "start", "end", "type", "info"]), file = out)
	for line in open(vcf, "r"):
		if line.startswith("#"):
			continue
		line = line.rstrip().split("\t")
		sv_type = line[7].split("SVTYPE=")[1].split(";")[0]

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

		print("\t".join([chr, str(st), str(fl), sv_type, line[7]]), file = out)

	out.close()


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
