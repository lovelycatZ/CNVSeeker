import sys 
import os


in_vcf = sys.argv[1]
out_vcf = sys.argv[2]
filtered_vcf = sys.argv[3]
data_type = sys.argv[4]



rd_tools = ["cn.mops", "CNVpytor", "gatk4", "ControlFREEC"]
non_rd_tools = ["delly", "Wham", "lumpy", "gridss", "Manta"]

with open(in_vcf, "r") as f, open(out_vcf, "w") as out, open(filtered_vcf, "w") as out1:
	tool_set = []
	for line in f:
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
		rd_supp_tool = []
		non_rd_supp_tool = []
		for tool_record in id_list.split(","):
			tool = tool_record.split("-")[0]
			if tool in rd_tools:
				rd_supp_tool.append(tool)
			elif tool in non_rd_tools:
				non_rd_supp_tool.append(tool)
			supp_tool.append(tool)

		if data_type == "LD_WGS":
			specifc_tools = ["lumpy", "gridss", "Wham"]
			large_supp_n = 3
		else:
			specifc_tools = ["lumpy"] if sv_type == "DEL" else ["Wham"]
			large_supp_n = 3

		if avg_len < 0:
			avg_len = -avg_len

		flag = 0
		if avg_len <= 50000:
			if (supp >= 3 and non_rd_supp_tool) or any(tool in supp_tool for tool in specifc_tools):
				flag = 1
		else:
			if len(rd_supp_tool) >= large_supp_n:
				flag = 1
				if all(x in rd_tools for x in supp_tool):
					st_list = []
					fl_list = []
					for tool_record in id_list.split(","):
						_, _, st_1, fl_1, _ = tool_record.strip().split("-")
						st_list.append(int(st_1))
						fl_list.append(int(fl_1))
					avg_start = min(st_list)
					avg_end = max(fl_list)
					avg_len = avg_end - avg_start + 1

			elif (len(non_rd_supp_tool) >= 3) or (len(rd_supp_tool) == 1 and len(non_rd_supp_tool) == 2) or (len(rd_supp_tool) == 2):
				if not os.path.exists("tmp/"):
					os.makedirs("tmp/")
				rd_tools1 = [x for x in rd_tools if x not in rd_supp_tool]
				threshold = 0.6
				flag1 = 0
				for tool in rd_tools1:
					os.system(f'echo -e "{chr}\t{avg_start}\t{avg_end}" > tmp/demo1.bed')
					os.system(f"grep '^#' {in_vcf} > tmp/demo2.vcf")
					os.system(f"grep {chr} {in_vcf} | grep '<{sv_type}>' | grep {tool} >> tmp/demo2.vcf")
					os.system(f"bedtools intersect -a tmp/demo1.bed -b tmp/demo2.vcf -wo > tmp/demo3.bed")
					with open("tmp/demo3.bed", "r") as f:
						total_length = 0
						id_list_new = [] 
						for line1 in f:
							if line1.startswith("#"):
								continue
							line1 = line1.strip().split("\t")
							supp = line1[10].split(";SUPP=")[-1].split(";")[0]
							id_list_1 = line1[10].split(";IDLIST=")[-1].split(";")[0]
							if supp == "1":
								id = line1[5]
								length = int(line1[-1])	
							else:
								for id1 in id_list_1.split(","):
									if id1.split("-")[0] == tool:
										id = id1
										break
								# st = max(avg_start, int(id.split("-")[2]))
								# fl = min(avg_end, int(id.split("-")[3]))
								st, fl = id.split("-")[2:4]
								length = int(fl) - int(st) + 1

							total_length += length
							id_list_new.append(id)
							id_list_new.append(f"{total_length}-{int(total_length) / int(avg_len)}")

					if int(total_length) / int(avg_len) >= threshold:
						id_list += "+" + ",".join(id_list_new)
						flag1 += 1
						if len(rd_supp_tool) == 2:
							if flag1 >= 1:
								flag = 1
								break
						else:
							if flag1 >= 2:
								flag = 1
								break
		
		line[1] = str(avg_start)
		if sv_type == "DEL":
			avg_len = -avg_len
		info = f"SVTYPE={sv_type};SVLEN={avg_len};END={avg_end};SUPP_TOOL={','.join(supp_tool)};ID_LIST={id_list}"
		line[7] = info
		if flag:
			print("\t".join(line), file = out)
		else:
			print("\t".join(line), file = out1)
