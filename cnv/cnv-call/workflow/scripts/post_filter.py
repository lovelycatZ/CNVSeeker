import pandas as pd
import shutil
import argparse
import os


################################
# Parse command-line arguments #
parser = argparse.ArgumentParser()
parser.add_argument("-e", "--exclude", type = str)
parser.add_argument("-i", "--input_vcf", type = str)
parser.add_argument("-o", "--output_vcf", type = str)
parser.add_argument("-p", "--overlap_threhold", type = float)
args = parser.parse_args()
################################


exclude = args.exclude 
input_vcf = args.input_vcf
output_vcf = args.output_vcf
overlap_threhold = args.overlap_threhold



def sort_vcf(vcf_file, output_file):
	header_lines = []
	variant_lines = []
	
	with open(vcf_file, 'r') as f:
		for line in f:
			if line.startswith('#'):
				header_lines.append(line)
			else:
				parts = line.strip().split('\t')
				chrom = parts[0]
				pos = int(parts[1])
				end = parts[7].split(";END=")[-1].split(";")[0]
				variant_lines.append((chrom, pos, end, line))
	
	variant_lines.sort(key=lambda x: (x[0], x[1], x[2]))
	
	with open(output_file, 'w') as f:
		f.writelines(header_lines)
		for _, _, _, line in variant_lines:
			f.write(line)



# 解析输入文件
header_lst = []
def vcf2bed(in_vcf, out_bed):
	with open(out_bed, "w") as out, open(in_vcf, "r") as vcf:
		for line in vcf:
			if line.startswith("#"):
				header_lst.append(line)
				continue

			line = line.rstrip().split("\t")
			chr = line[0]
			st = line[1]
			info = line[7]

			if not chr.startswith("chr"):
				chr = "chr" + chr
			
			if info.startswith("END="):
				fl = info.split(";")[0].split("=")[-1]
			elif ";END=" in info:
				fl = info.split(";END=")[1].split(";")[0]
			elif "SVLEN=" in info:
				sv_len = int(info.split("SVLEN=")[1].split(";")[0])
				if sv_len < 0:
					sv_len = -sv_len
				fl = str(int(st) + sv_len - 1)
			else:
				continue
			
			cnv_len = str(int(fl) - int(st) + 1)
			cnv = "\t".join([chr, st, fl] + line + [cnv_len])
			print(cnv, file = out)


tmp_dir = "./tmp/"
if not os.path.exists(tmp_dir):
	os.makedirs(tmp_dir)
tmp_bed = tmp_dir + "tmp.bed"
tmp_out_vcf = tmp_dir + "tmp_out.vcf"
intersections = tmp_dir + "intersections.bed"
no_intersections = tmp_dir + "no_intersections.bed"
merged_intersections = tmp_dir + "merged_intersections.bed"


vcf2bed(input_vcf, tmp_bed)
cmd1 = f"bedtools intersect -a {tmp_bed} -wo -b {exclude} > {intersections}"
cmd2 = f"bedtools intersect -a {tmp_bed} -v -b {exclude} > {no_intersections}"
os.system(cmd1)
os.system(cmd2)


if os.path.getsize(intersections):
	# 读入数据框并按CNV分组
	df = pd.read_csv(intersections, sep = "\t", low_memory = False, header = None)
	num_columns = df.shape[1]

	with open(exclude, 'r') as f:
		exclude_col_len = len(f.readline().strip().split("\t"))

	merged_df = df.groupby([i for i in range(3, num_columns - exclude_col_len - 1)])[num_columns - 1].sum().reset_index()
	merged_df.to_csv(merged_intersections, sep = "\t", index = False, header = False)

	with open(tmp_out_vcf, "w") as out_vcf:
		# 写入 vcf header
		for line in header_lst:
			out_vcf.write(line)

		# 过滤与 exclude 重叠超过阈值的 CNV 并写入文件
		for _, row in merged_df.iterrows():
			row = list(row)
			cnv_len = int(row[-2])
			overlap_len = int(row[-1])
			if overlap_len / cnv_len <= float(overlap_threhold):
				out_lst = [str(item) for item in row[0:-2]]
				print("\t".join(out_lst), file = out_vcf)

		# 写入与 exclude 无重叠的CNV
		for line in open(no_intersections):
			line = line.rstrip().split("\t")[3:-1]
			print("\t".join(line), file = out_vcf)

	# 对 VCF 文件排序
	sort_vcf(tmp_out_vcf, output_vcf)
	
else:
	os.system(f"cp {input_vcf} {output_vcf}")

# shutil.rmtree(tmp_dir)
