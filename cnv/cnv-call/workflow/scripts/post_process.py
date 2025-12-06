def cnmops_postprocess(input_file, sample_name, output_file):
	with open(input_file, 'r') as input_tsv, open(output_file, 'w') as vcf:
		vcf.write('##fileformat=VCFv4.2\n')
		vcf.write('##source=cn.MOPS\n')
		vcf.write('##ALT=<ID=DEL,Description="Deletion">\n')
		vcf.write('##ALT=<ID=DUP,Description="Duplication">\n')
		vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
		vcf.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
		vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n')
		vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
		vcf.write('##FORMAT=<ID=CN,Number=1,Type=String,Description="Copy number value">\n')
		vcf.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name + "_cn.MOPS"}\n')
		for line in input_tsv:
			out_lst= []
			if line.startswith("seqnames"):
				continue
			
			line = line.strip().split('\t')
			if sample_name not in line[6]:
				continue

			out_lst.append(line[1])
			out_lst.append(line[2])
			out_lst.append(".")
			out_lst.append("N")
			cn = int(line[-1][-1])
			if cn < 2:
				cnv_type = "<DEL>"
				if cn == 0:
					gt = "1/1"
				else:
					gt = "0/1"
			else:
				cnv_type = "<DUP>"
				if cn == 3:
					gt = "0/1"
				else:
					# we can't really say if this is a hom dup or het dup with higher copy number.
					gt = "./1"
			out_lst.append(cnv_type)
			out_lst.append(".")
			out_lst.append("PASS")
			info = f"SVTYPE={cnv_type[1:-1]};SVLEN={line[4]};END={line[3]}"
			out_lst.append(info)
			out_lst.append("GT:CN")
			out_lst.append(f"{gt}:{cn}")
			outline = "\t".join(out_lst)
			print(outline, file = vcf)


# def CNVpytor_postprocess(input_file, sample_name, output_file, p1_thd, q0_thd):
# 	with open(input_file, "r") as i, open(output_file, "w") as o:
# 		for line in i:		
# 			if line.startswith("##"):
# 				o.write(line)
# 				continue
# 			if line.startswith("#CHROM"):
# 				line = line.rstrip().split("\t")[0:9]
# 				line.append(sample_name + "_CNVpytor")
# 				print("\t".join(line), file = o)
# 				continue
# 			line = line.rstrip().split("\t")
# 			p1 = line[7].split("pytorP1=")[-1].split(";")[0]
# 			if p1 == "nan":
# 				continue
# 			# p1 = eval(p1)
# 			# # p2 = eval(line[7].split("natorP2=")[-1].split(";")[0])
# 			# q0 = eval(line[7].split("pytorQ0=")[-1].split(";")[0])
# 			# if p1 >= eval(p1_thd) or q0 >= eval(q0_thd):
# 			# 	continue
# 			gt = line[9].split(":")[0]
# 			line[7] += f";GT={gt}"
# 			print("\t".join(line), file = o)


# def delly_postprocess(input_file, sample_name, output_file):
# 	with open(input_file, 'r') as input_vcf, open(output_file, 'w') as vcf:
# 		for line in input_vcf:
# 			if line.startswith("##"):
# 				vcf.write(line)
# 				continue
# 			if line.startswith("#CHROM"):
# 				line = line.strip().split('\t')[0:9]
# 				line.append(sample_name + "_delly")
# 				print("\t".join(line), file = vcf)
# 				continue
# 			line = line.strip().split('\t')
# 			if line[6] != "PASS":
# 				continue
# 			gt = line[9].split(":")[0]
# 			line[7] += f";GT={gt}"
# 			print("\t".join(line), file = vcf)


def ExomeDepth_postprocess(input_file, sample_name, output_file):
	with open(input_file, 'r') as input_vcf, open(output_file, 'w') as vcf:
		vcf.write('##fileformat=VCFv4.2\n')
		vcf.write('##source=ExomeDepth\n')
		vcf.write('##ALT=<ID=DEL,Description="Deletion">\n')
		vcf.write('##ALT=<ID=DUP,Description="Duplication">\n')
		vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
		vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
		vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n')
		vcf.write('##INFO=<ID=BF,Number=1,Type=Float,Description="Bayes Factor">\n')
		vcf.write('##INFO=<ID=EXONS,Number=.,Type=Integer,Description="Number of exons affected">\n')
		vcf.write('##INFO=<ID=READSEXP,Number=.,Type=Integer,Description="Number of reads expected">\n')
		vcf.write('##INFO=<ID=READSOBS,Number=.,Type=Integer,Description="Number of reads observed">\n')
		vcf.write('##INFO=<ID=READSRATIO,Number=.,Type=Integer,Description="Ratio expected vs. observed">\n')
		vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
		vcf.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name + "_ExomeDepth"}\n')
		for line in input_vcf:
			if "start.p" in line:
				continue
			line = line.strip().split('\t')
			chr = line[6]
			pos = line[4]
			end = line[5]
			id = "."
			ref = "N"
			alt = '<DEL>' if line[2] == "deletion" else '<DUP>'
			qual = "."
			filter = "PASS"
			if eval(line[11]) > 0.5 and eval(line[11]) < 1.5:
				gt = "0/1"  
			else: 
				gt = "1/1"
			type = alt[1:-1]
			len = int(end) - int(pos) + 1
			info = f'SVTYPE={type};END={end};SVLEN={len};EXONS={line[3]};READSEXP={line[9]};READSOBS={line[10]};READSRATIO={line[11]};BF={line[8]};GT={gt}'
			outline = f'{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{"GT"}\t{gt}'
			print(outline, file = vcf)
			

def gatk4_postprocess(input_file, sample_name, output_file):
	with open(input_file, 'r') as input_vcf, open(output_file, 'w') as vcf:
		for line in input_vcf:
			if line.startswith("##"):
				vcf.write(line)
				continue

			if line.startswith("#CHROM"):
				line = line.rstrip().split('\t')[0:9]
				line.append(sample_name + "_gatk4")
				print("\t".join(line), file = vcf)
				continue

			lst = line.strip().split('\t')
			if lst[4] == ".":
				continue

			vcf.write(line)
			# gt = line[9][0]
			# cn = line[9][2]
			# if gt == "1":
			# 	line[4] = '<DEL>'
			# 	sv_type = "DEL"
			# elif gt == "2":
			# 	line[4] = '<DUP>'
			# 	sv_type = "DUP"
			# else:
			# 	continue			

			# if gt == "1" and cn == 0:
			# 	GT = "1/1"
			# elif gt == "1" and cn == 1:
			# 	GT = "0/1"
			# elif gt == "2" and cn == 3:
			# 	GT = "0/1"
			# else:
			# 	GT = "./1"
			
			# line[7] += f";SVTYPE={sv_type};GT={GT}"
			# print("\t".join(line), file = vcf)
			
			

# def lumpy_postprocess(input_file, sample_name, output_file):
# 	with open(input_file, 'r') as input_vcf, open(output_file, 'w') as vcf:
# 		for line in input_vcf:
# 			if line.startswith("##"):
# 				vcf.write(line)
# 				continue
# 			if line.startswith("#CHROM"):
# 				line = line.rstrip().split('\t')[0:9]
# 				line.append(sample_name + "_lumpy")
# 				print("\t".join(line), file = vcf)
# 				continue
# 			line = line.rstrip().split('\t')
# 			alt = line[4]
# 			if alt not in ["<DUP>", "<DEL>"]:
# 				continue
# 			gt = line[9].split(":")[0]
# 			line[7] += f";GT={gt}"
# 			print("\t".join(line), file = vcf)

			

def xhmm_postprocess(input_file, sample_name, output_file):
	with open(input_file, 'r') as input_vcf, open(output_file, 'w') as vcf:
		for line in input_vcf:
			if line.startswith("##"):
				vcf.write(line)
				continue

			line = line.strip().split('\t')
			if "#CHROM" in line:
				if sample_name in line:
					col_index = line.index(sample_name)
				else:
					print(f"{sample_name} had been filtered out by xhmm")
					break
				line = line[0:9]
				line.append(sample_name + "_xhmm")
				print("\t".join(line), file = vcf)
				continue
			
			info = line[col_index]
			if info[-1] == "Y":
				line = line[0:9]
				line.append(info)
				gt = info[0]
				# rd = info.split(":")[-3]
				if gt == "1":
					alt = '<DEL>'
				elif gt == "2":
					alt = '<DUP>'
				else:
					continue
				
				line[4] = alt
	
				# GT = ""
				# if gt == "1" and float(rd) < 0.25:
				# 	GT = "1/1"
				# elif gt == "1" and float(rd) >= 0.25:
				# 	GT = "0/1"
				# elif gt == "2" and float(rd) <= 1.75:
				# 	GT = "0/1"
				# elif gt == "2" and float(rd) > 1.75:
				# 	GT = f"./1"
				# else:
				# 	GT = "./.:."
				
				# line[7] += f";GT={GT}"
				print("\t".join(line), file = vcf)
