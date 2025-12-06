import sys 
import os


tmp_vcf = sys.argv[1]
out_vcf = sys.argv[2]
filtered_vcf = sys.argv[3]
data_type = sys.argv[4]
exons = sys.argv[5]


base_path = os.path.dirname(tmp_vcf)
base_name = os.path.basename(tmp_vcf).split(".")[0]
tmp_dir = base_path + f"/tmp_{base_name}"
if not os.path.exists(tmp_dir):
	os.makedirs(tmp_dir)

with open(tmp_vcf, "r") as f, open(out_vcf, "w") as out, open(filtered_vcf, "w") as out1:
	tool_set = []
	for i, line in enumerate(f):
		if line.startswith("#"):
			out.write(line)
			continue
		
		line = line.strip().split("\t")
		chr = line[0]
		info = line[7]
		avg_start = int(info.split("AVG_START=")[-1].split(";")[0].split(".")[0])
		avg_end = int(info.split("AVG_END=")[-1].split(";")[0].split(".")[0])
		avg_len = avg_end - avg_start + 1
		sv_type = info.split("SVTYPE=")[-1].split(";")[0]
		supp = int(info.split(";SUPP=")[-1].split(";")[0])
		id_list = info.split(";IDLIST=")[-1].split(";")[0]

		supp_tool = []
		for tool_record in id_list.split(","):
			tool = tool_record.split("-")[0]
			supp_tool.append(tool)

		# 统计外显子个数
		tmp_CNV = tmp_dir + f"/tmp.cnv.{i}.bed"
		tmp_CNV_exons = tmp_dir + f"/tmp.cnv.exons.{i}.bed"
		with open(tmp_CNV, "w") as f:
			print("\t".join([chr, str(avg_start), str(avg_end)]), file = f)
		os.system(f"bedtools intersect -a {tmp_CNV} -b {exons} > {tmp_CNV_exons}")
		with open(tmp_CNV_exons, 'r') as file:
			lines = file.readlines()
		exons_count = len(lines)

		if avg_len < 0:
			avg_len = -avg_len

		flag = 0
		if sv_type == "DEL":
			if all(tool in supp_tool for tool in ["ExomeDepth", "cn.mops"]) or any(tool in supp_tool for tool in ["gatk4", "ECOLE"]):
				flag = 1
			
		else:
			if all(tool in supp_tool for tool in ["ExomeDepth", "cn.mops", "xhmm"]) or "ECOLE" in supp_tool:
				flag = 1
				   
		if not flag:
			continue

		line[1] = str(avg_start)
		info = f"SVTYPE={sv_type};SVLEN={avg_len};END={avg_end};SUPP_TOOL={','.join(supp_tool)};ID_LIST={id_list};EXONS={exons_count}"
		line[7] = info
		if flag:
			print("\t".join(line), file = out)
		else:
			print("\t".join(line), file = out1)

#os.system(f"rm -rf {tmp_dir}")