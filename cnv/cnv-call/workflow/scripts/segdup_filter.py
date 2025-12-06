import argparse
import pybedtools as pbt
import pandas as pd
import os


################################
# Parse command-line arguments #
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_vcf", type = str)
parser.add_argument("-e", "--exclude", type = str)
parser.add_argument("-o", "--output_vcf", type = str)
# parser.add_argument("-n", "--name", type = str)
parser.add_argument("-c", "--column_nums", type = int)
parser.add_argument("-p", "--overlap_threhold", type = float)
args = parser.parse_args()
################################


exclude = args.exclude 
input_vcf = args.input_vcf
output_vcf = args.output_vcf
# data_name = args.name
col_num = args.column_nums
overlap_threhold = args.overlap_threhold


# 提取 vcf 文件以 #CHROM 开头的行
vcf_header = ""
with open(input_vcf, "r") as f:
	for line in f:
		if line.startswith("#"):
			vcf_header += line
		if line.startswith("#CHROM"):
			header_line = line.strip()
			break


vcf_header_lst = header_line.split("\t")
bed_header_lst = ["CHR", "START", "END"]
header_lst = vcf_header_lst + bed_header_lst + ["overlap_len"]

df_exclude = pd.read_csv(exclude, sep = "\t", header = None, comment = "#", usecols = [0,1,2], names = bed_header_lst)
bedtool_ex = pbt.BedTool.from_dataframe(df_exclude)
bedtool_vcf = pbt.BedTool(input_vcf)

inter = bedtool_vcf.intersect(bedtool_ex, wao = True)
df_inter = inter.to_dataframe(names = header_lst)

df_inter.drop(columns = bed_header_lst, inplace = True)
# 如果overlap_len列是 NaN, 填充为0
df_inter['overlap_len'] = df_inter['overlap_len'].fillna(0)


if df_inter.empty:
	os.system(f"cp {input_vcf} {output_vcf}")
else:
	df_inter = df_inter[df_inter["ALT"].str.contains("DUP") | df_inter["ALT"].str.contains("DEL")].copy()
	df_inter = df_inter.groupby(vcf_header_lst)["overlap_len"].sum().reset_index()
	df_inter['length'] = df_inter['INFO'].apply(lambda x: abs(int(x.split("SVLEN=")[1].split(";")[0])))
	df_inter['overlap'] = df_inter["overlap_len"] / df_inter["length"]
	with open(output_vcf, "w") as out:
		out.write(vcf_header)
		df_inter = df_inter[(df_inter["ALT"].str.contains("DEL")) | (df_inter["overlap"] <= overlap_threhold) | (df_inter['length'] < 5000)].copy()
		df_inter.drop(columns = ["length", "overlap", "overlap_len"], inplace = True)
		df_inter.to_csv(out, sep = "\t", header = False, index = False)
