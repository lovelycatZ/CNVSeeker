import os
import sys


input_vcf = sys.argv[1]
output_vcf = sys.argv[2]
tool = sys.argv[3]
sample_id = sys.argv[4] 
remain_sv_type = sys.argv[5]  # "DEL" or "DUP" or "ALL"
min_len = int(sys.argv[6])
max_len = int(sys.argv[7])
apply_filter = sys.argv[8]
sex = sys.argv[9]


def pass_filter(tool, record_lst):

	gt = record_lst[-1].split(":")[0]
	if gt == "0/0":
		return False
	
	info = record_lst[7]
	if tool == "CNVpytor":
		pytorRD = info.split(";pytorRD=")[-1].split(";")[0]
		PN = info.split(";pN=")[-1].split(";")[0]
		DG = info.split(";dG=")[-1].split(";")[0]
		P1 = info.split(";pytorP1=")[-1].split(";")[0]
		Q0 = info.split(";pytorQ0=")[-1].split(";")[0]
		if P1 == "nan":
			return False
		if eval(P1) > 0.05:
			return False
		# if float(PN) > 0.6:
		# 	return False
		# if int(DG) <= 100000:
		# 	return False
		if eval(Q0) > 0.9:
			return False
		
		if eval(pytorRD) > 0.75 and eval(pytorRD) < 1.25:
			print("likely copy number neutral, no CNV")
			return False

	elif tool == 'lumpy':
		qual = eval(record_lst[5])
		if qual < 120:
			return False

	elif tool == "delly":
		if "PASS" not in record_lst[6]:
			return False
		if "PASS" not in record_lst[9]:
			return False
		if "0/0" in record_lst[9]:
			return False
	
	elif tool == "Manta":
		if "PASS" not in record_lst[6]:
			return False	
		if "PASS" not in record_lst[9]:
			return False
			
	elif tool == "Wham":
		sv_type = record_lst[4][1:-1]
		TF = info.split(";TF=")[-1].split(";")[0]
		U = info.split(";U=")[-1].split(";")[0]
		CW_list = info.split(";CW=")[-1].split(";")[0].split(",")
		CW_list = [eval(i) for i in CW_list]
		max_CW = max(CW_list)
		if sv_type == "DEL" and  eval(TF) < 3:
			return False
		if sv_type == "DUP" and eval(U) < 3:
			return False
		if max_CW < 0.2:
			return False

	elif tool == "ExomeDepth":
		svtype = record_lst[4][1:-1]
		BF = info.split(";BF=")[-1].split(";")[0]
		ratio = info.split(";READSRATIO=")[-1].split(";")[0]
		if eval(BF) < 15:
			return False
		if (svtype == "DUP" and eval(ratio) < 1.1) or (svtype == "DEL" and eval(ratio) > 0.8):
			return False
	
	elif tool == "gatk4":
		qual = record_lst[5]
		gt_info = record_lst[-1]
		cn = eval(gt_info.split(":")[1])
		nt = eval(gt_info.split(":")[2])

		if cn == 0:
			threshold = min(1000, max(400, 10*nt))
		elif cn == 1:
			threshold = min(1000, max(100, 10*nt))
		else:
			threshold = min(400, max(50, 4*nt))
			# threshold = 1000

		if eval(qual) < threshold:
			return False

	elif tool == "ControlFREEC":
		W = info.split(";W=")[-1].split(";")[0]
		KS = info.split(";KS=")[-1].split(";")[0]
		if eval(W) > 0.05 or eval(KS) > 0.05:
			return False

	# TGS
	elif tool in ["cuteSV", "SVision_pro", "Sniffles", "DeBreak", "NanoVar"]:
		if "PASS" not in record_lst[6]:
			return False
		if any(x in record_lst[9] for x in ["0/0", "./."]):
			return False
		
	else:
		pass

	return True


def standardize(input_vcf, output_vcf, remain_sv_type, tool, sample_id, min_len, max_len, apply_filter, sex):
	CHROM_LIST = ["chr" + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
	remain_sv_type_list = remain_sv_type.split(',')
	if sex in ['XY', 'M', 'Male', 'male', 'MALE', 'm']:
		sex = 'M'
	else:
		sex = 'F'

	with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
		for line in infile:
			# 跳过空行
			if not line.strip():
				continue

			if line.startswith('#'):
				if "##contig=" in line:
					chr = line.split('ID=')[1].split(',')[0]
					if chr not in CHROM_LIST:
						continue
				if line.startswith("#CHROM"):
					header_parts = line.strip().split('\t')
					header_parts[-1] = sample_id
					line = '\t'.join(header_parts) + '\n'
				outfile.write(line)

			else:
				columns = line.strip().split('\t')

				chrom = columns[0]
				if chrom not in CHROM_LIST:
					continue
				
				if sex == 'F' and chrom == 'chrY':
					continue
				
				alt = columns[4]
				info_field = columns[7]
				info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)

				sv_type = info_dict.get("SVTYPE", "")
	
				if sv_type not in remain_sv_type_list:
					# 尝试从ALT字段中获取SVTYPE
					if "DEL" in alt:
						sv_type = "DEL"
					elif "DUP" in alt:
						sv_type = "DUP"
					else:
						continue
				
				if sv_type not in remain_sv_type_list:
					continue
				
				# 修复 xhmm 的 SVTYPE
				if tool == 'xhmm':
					info_dict['SVTYPE'] = sv_type

				if apply_filter:
					if not pass_filter(tool, columns):
						continue

				start = int(columns[1])
				if start == 0:
					start = 1
					columns[1] = str(start)

				if "SVLEN" not in info_dict:
					end = info_dict.get("END", None)
					end = int(end)
					sv_len = end - start + 1
					if sv_type == 'DEL':
						sv_len = -sv_len
					info_dict['SVLEN'] = str(sv_len)

				if "END" not in info_dict:
					sv_len = info_dict.get('SVLEN', None)
					end = start + abs(int(sv_len)) - 1
					info_dict['END'] = str(end)
				
				if "SVTYPE" not in info_dict:
					info_dict['SVTYPE'] = sv_type

				sv_len = abs(int(info_dict.get('SVLEN', None)))
				if sv_len < min_len or sv_len > max_len:
					continue
				
				# 重新生成ID字段
				end = info_dict.get("END")
				id = f"{sample_id}-{tool}-{chrom}-{start}-{end}-{sv_type}"
				columns[2] = id

				# 把info_dict重新组装回info字段
				new_info_field = ';'.join([f"{key}={value}" for key, value in info_dict.items()])
				columns[7] = new_info_field

				print("\t".join(columns), file=outfile)


if __name__ == "__main__":
	standardize(input_vcf, output_vcf, remain_sv_type, tool, sample_id, min_len, max_len, apply_filter, sex)
