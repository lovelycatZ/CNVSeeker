import time
t1 = time.perf_counter()

from interpretation.resources import *
from interpretation.annotate import *
from interpretation.utils import *
from tqdm import tqdm
import copy
import os

t2 = time.perf_counter()
print('Elapsed time of importing packages:', '{0:.1f}'.format(t2-t1), 'seconds')


def annotation(CNV):
	cnv = CNVinfo(CNV[0], CNV[1], CNV[2], CNV[3], CNV[-1])
	rules = {} # 记录该CNV所满足的证据项	
	if cnv.type == "DEL":
		details = copy.deepcopy(DETAIL_TABLE["DEL"])
	else:
		details = copy.deepcopy(DETAIL_TABLE["DUP"])

	# 选择参考基因组版本
	if cnv.build == "hg19":
		DB_gene_info = DB_gene_info_hg19
		DB_cytoband = DB_cytoband_hg19
		DB_ClinGen_regions = DB_ClinGen_regions_hg19
		DB_ClinVar = DB_ClinVar_hg19
		DB_Clinvar_null_variant = DB_ClinVar_null_variant_hg19
		DB_gnomAD_SV = DB_gnomAD_SV_hg19
		DB_DGV_GS_CNV = DB_DGV_GS_CNV_hg19
		DB_DDD_CNV = DB_DDD_CNV_hg19
		DB_ClinVar_CNV = DB_ClinVar_CNV_hg19
		#DB_dbvar_common = DB_dbvar_common_hg19
		#DB_dbvar_pathogenic = DB_dbvar_pathogenic_hg19
	elif cnv.build == "hg38":
		DB_gene_info = DB_gene_info_hg38
		DB_cytoband = DB_cytoband_hg38
		DB_ClinGen_regions = DB_ClinGen_regions_hg38
		DB_ClinVar = DB_ClinVar_hg38
		DB_Clinvar_null_variant = DB_ClinVar_null_variant_hg38
		DB_gnomAD_SV = DB_gnomAD_SV_hg38
		DB_DGV_GS_CNV = DB_DGV_GS_CNV_hg38
		DB_DDD_CNV = DB_DDD_CNV_hg38
		DB_ClinVar_CNV = DB_ClinVar_CNV_hg38
		#DB_dbvar_common = DB_dbvar_common_hg38
		#DB_dbvar_pathogenic = DB_dbvar_pathogenic_hg38

	out_record = OutRecord()
	#out_record.chr = cnv.chr
	#out_record.start = cnv.start
	#out_record.end = cnv.end
	#out_record.type = cnv.type
	out_record.build = cnv.build
	out_record.size = cnv.size()
	out_record.cytoband = CNV_cytoband(cnv, DB_cytoband)

	gene_dict = count_overlap_genes(cnv, DB_gene_info)
	procoding_genes_nums = len(gene_dict)
	gene_str = ", ".join(sorted(gene_dict.keys()))
	
	out_record.gene_count = procoding_genes_nums

	## 第1部分评分
	if procoding_genes_nums != 0:
		rules["1A"] = True
		if procoding_genes_nums > 1:
			details["1A"] += f"This CNV overlaps with {procoding_genes_nums} protein-coding genes including: \n {gene_str}."
		else:
			details["1A"] += f"This CNV fully contained within the protein-coding gene: \n {gene_str}."	
	else:
		rules["1B"] = True
		details["1B"] += "This CNV DOES NOT contain any protein-coding genes or known functionally important elements."

	## 评估 ClinGen 的致病性/良性基因
	# 记录 CNV 是否与其他基因有部分重叠
	cnv_gene_partial_overlaps_count = 0
	gene_list = []
	result = DB_gene_info.overlap(cnv.chr, cnv.start, cnv.end)
	for record, overlap, coverage in result:
		row = record._asdict()
		gene_chr = row["gene_chr"]
		gene_st = int(row["gene_start"]) + 1
		gene_fl = row["gene_end"]
		overlap_mode = check_overlap(overlap, coverage)
		if overlap_mode == "1":
			cnv_gene_partial_overlaps_count += 1 

		gene_symbol = row["gene_symbol"]
		out_record.genes.append(f"{gene_symbol}:{gene_chr}:{gene_st}-{gene_fl}({overlap_mode})")

		gene_group_id = row["gene_group_id"]
		HI_score = row["HI_Score"]
		TS_score = row["TS_Score"]
		morbid_phenotype = row["morbid_phenotype"]

		gene_record = [gene_symbol, [], ""]
		if gene_group_id != "-":
			gene_record[1] = str(gene_group_id).split("|")

		if morbid_phenotype != "-":
			gene_record[2] = morbid_phenotype
		gene_list.append(gene_record)
		
		## 判断基因的病态性
		if morbid_phenotype != "-":			
			out_record.OMIM_genes.append(f"{gene_symbol}:{gene_group_id}:[{morbid_phenotype}]")
		
		## 判断基因的剂量敏感性
		if HI_score not in ["Not yet evaluated", "-"]:
			HI_score = int(HI_score)
			# if HI_score in [3, 2, 40]:
			out_record.HI_genes.append(f"{gene_symbol}:HI-{HI_score}")

		if TS_score not in ["Not yet evaluated", "-"]:
			TS_score = int(TS_score)
			# if TS_score in [3, 2, 40]:				
			out_record.TS_genes.append(f"{gene_symbol}:TS-{TS_score}")

		HI_index = row["HI_index"]
		pLI_score = row["pLI"]
		oe_lof_upper = row["oe_lof_upper"]
		out_record.pLI_loeuf_HI.append(f"{gene_symbol}(pLI={pLI_score};loeuf={oe_lof_upper};HI%={HI_index})")
		
		if all(x != "-" for x in [HI_index, pLI_score, oe_lof_upper]):
			if eval(HI_index) < 10 and eval(pLI_score) > 0.9 and eval(oe_lof_upper) < 0.35:
				out_record.p_HI_genes.append(gene_symbol)

				## DEL - 2H
				if "DEL" in cnv.type and HI_score != 3:
					rules["2H"] = True
					if "2H" not in details:
						details["2H"] = f"Genes overlapping with the CNV that are predicted to be haploinsufficient: \n"
					
					details["2H"] += f"""* {gene_symbol}: ClinGen Haploinsufficiency score: {HI_score}. \n gnomAD PLI score: {pLI_score} \
(the upper bound of the observed/expected confidence interval: {oe_lof_upper}). \n DECIPHER HI score: {HI_index}.
"""				

		# if all(score not in [3, 2, 40] for score in [HI_score, TS_score]):
		# 	continue
		
		p_LOF_threshold = 500
		null_variant_nums = gene_null_variant([gene_chr, gene_st, gene_fl], DB_Clinvar_null_variant)
		if null_variant_nums >= 500:
			out_record.p_LOF_genes.append(f"{gene_symbol}(null variants: {null_variant_nums})")
			if HI_score not in [3, 2]:
				HI_score = 3
		elif pLI_score != "-" and null_variant_nums >= 200:
			if eval(pLI_score) < 0.9:
				continue
			out_record.p_LOF_genes.append(f"{gene_symbol}(null variants: {null_variant_nums})")
			if HI_score not in [3, 2]:
				HI_score = 3


		#### DUP ####
		## 与 TS 基因有交集
		if "DUP" in cnv.type:
			if TS_score in [3, 2]:
				if overlap_mode in [">", "="]:
					rules[f"2A_TS-{TS_score}"] = True
					if "* Genes: \n" not in details["2A"]:
						details["2A"] += "* Genes: \n"
					details["2A"] += f" - {gene_symbol}: ClinGen Triplosensitivity score: {TS_score}. \n"
					"""
					需要增加HI基因的信息
					"""

			elif TS_score == 40:
				if overlap_mode in ["<", "="]:
					rules["2C"] = True
					if "* Genes: \n" not in details["2C"]:
						details["2C"] += "* Genes: \n"
					details["2C"] += f" - {gene_symbol}: ClinGen Triplosensitivity score: {TS_score}. \n "

			## 与 HI 基因有交集
			if HI_score in [3, 2] and overlap_mode == "<":
				if overlap_mode in "<":
					tx_refseq_name = row["RefSeq_tx_id"].split("|")[0]
					pvs1_strength, pvs1_class = PVS1_Strength(cnv, tx_refseq_name)
					rules[f'2I-{pvs1_class}_HI-{HI_score}'] = True
					details[f'2I-{pvs1_class}_HI-{HI_score}'] = f""" PVS1 strength is {pvs1_strength} according to ClinGen SVI working group PVS1 specifications.
* ClinGen Haploinsufficiency score: {HI_score}.
* The MANE Select transcript ID: {tx_refseq_name}.
* gnomAD PLI score: {pLI_score}.
* the upper bound of the observed/expected confidence interval: {oe_lof_upper}).
* DECIPHER HI score: {HI_index}.
"""

				elif overlap_mode == "1":
					## should be associated with patient’s phenotype
					patient_phenotype = 0
					## patient’s phenotype is either inconsistent with what is expected for LOF of that gene OR unknown.
					if patient_phenotype == 0:
						rules['2J'] = True
						details["2J"] += f"* {gene_symbol}: ClinGen Triplosensitivity score: {HI_score}. \n"
					
					else:
					## patient’s phenotype is highly specific and consistent with what is expected for LOF of that gene.
						rules['2K'] = True
						details["2K"] += f"* {gene_symbol}: ClinGen Triplosensitivity score: {HI_score}. \n"

			# 2L
			if (HI_score not in [3, 2]) and (TS_score not in [3, 2]):
				rules["2L"] = True
				details["2L"] += f"* {gene_symbol}: ClinGen Haploinsufficiency score: {HI_score}. \n ClinGen Triplosensitivity score: {HI_score}. \n "
		#### DUP ####	
	

		#### DEL ####
		else:
			# if HI_score in [3, 2, 30]:
			if HI_score in [3, 2]:
				## 2A: 和已知HI基因完全重叠
				if overlap_mode in [">", "="]:
					rules[f'2A_HI-{HI_score}'] = True
					if "* Genes: \n" not in details["2A"]:
						details["2A"] += "* Genes: \n"					
					details["2A"] += f" - {gene_symbol}: ClinGen Haploinsufficiency score: {HI_score}.\n "
					"""
					需要增加HI基因的信息
					"""
					
				## 和已知HI基因部分重叠
				elif overlap_mode == "1":
					HI_gene_info = gene_content(row, DB_ClinVar)
					strand, five_prime_UTR, three_prime_UTR, cds, exons, snvs_str = HI_gene_info

					## 基因位于+链
					if strand == "+":
						## 2C
						if cnv.end <= gene_fl:
							## 2C-1
							if cnv.end >= cds[0]:
								rules[f'2C-1_HI-{HI_score}'] = True
								details["2C-1"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV spans the 5' UTR and coding sequence of {gene_symbol}.\n "
							
							## 2C-2
							if cnv.end >= five_prime_UTR[0] and cnv.end <= five_prime_UTR[1]:
								rules[f'2C-2_HI-{HI_score}'] = True
								details[f"2C-2"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV overlaps with the 5' UTR of {gene_symbol} only.\n "
								
						## 2D
						else:
							## 2D-1
							if cnv.start >= three_prime_UTR[0] and cnv.start <= three_prime_UTR[1]:
								rules[f'2D-1_HI-{HI_score}'] = True
								details["2D-1"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV overlaps the 3' UTR of {gene_symbol} only.\n "

							if len(exons) == 1:
								if cnv.start <= exons[-1][0]:
									if snvs_str:
										# snv_count = 0
										# for snv in snv_list:
										# 	if int(snv[1]) >= cnv.start:
										# 		snv_count += 1
										# ## 2D-2
										# if snv_count:
										rules[f'2D-2_HI-{HI_score}'] = True
										details["2D-2"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV overlaps the only (last) exon of {gene_symbol} and {snvs_str} established pathogenic small variants have been report in this exon incluing {snvs_str}. \n "

									else:
										## 2D-3
										rules[f'2D-3_HI-{HI_score}'] = True
										details["2D-3"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV overlaps the only (last) exon of {gene_symbol} and No established pathogenic small variants have been report in this exon. \n "

							else:								
								if cnv.start >= exons[-2][1] and cnv.start <= exons[-1][0]:
									if snvs_str:
										# snv_count = 0
										# for snv in snv_list:
										# 	if int(snv[1]) >= cnv.start:
										# 		snv_count += 1
										# ## 2D-2
										# if snv_count:
										rules[f'2D-2_HI-{HI_score}'] = True
										details["2D-2"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV overlaps the last exon of {gene_symbol} only and {snvs_str} established pathogenic small variants have been report in this exon incluing {snvs_str}. \n "

									else:
										## 2D-3
										rules[f'2D-3_HI-{HI_score}'] = True
										details["2D-3"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV overlaps the last exon of {gene_symbol} only and No established pathogenic small variants have been report in this exon. \n "

								if cnv.start <= exons[-2][0]:
									## 2D-4
									rules[f'2D-4_HI-{HI_score}'] = True
									details["2D-4"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: plus (+) strand. \n This CNV spans the 3' UTR and includes other exons in addition to the last exon. Nonsense-mediated decay expected to occur. \n "

					## 基因位于-链
					else:
						## 2C
						if cnv.start >= gene_st:
							## 2C-1
							if cnv.start <= cds[-1]:
								rules[f'2C-1_HI-{HI_score}'] = True
								details["2C-1"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV spans the 5' UTR and coding sequence of {gene_symbol}.\n "

							## 2C-2
							if cnv.start >= five_prime_UTR[0] and cnv.start <= five_prime_UTR[1]:
								rules[f'2C-2_HI-{HI_score}'] = True
								details["2C-2"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV overlaps with the 5' UTR of {gene_symbol} only.\n "
								
						## 2D
						else:
							## 2D-1
							if cnv.end >= three_prime_UTR[0] and cnv.end <= three_prime_UTR[1]:
								rules[f'2D-1_HI-{HI_score}'] = True
								details["2D-1"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV overlaps the 3' UTR of {gene_symbol} only.\n "

							if len(exons) == 1:
								if cnv.end >= exons[0][1]:
									if snvs_str:
										# snv_count = 0
										# for snv in snv_list:
										# 	if int(snv[1]) <= cnv.end:
										# 		snv_count += 1
										# ## 2D-2
										# if snv_count:
										rules[f'2D-2_HI-{HI_score}'] = True
										details["2D-2"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV overlaps the only (last) exon of {gene_symbol} and {snvs_str} established pathogenic small variants have been report in this exon incluing {snvs_str}. \n "

									else:
										## 2D-3
										rules[f'2D-3_HI-{HI_score}'] = True
										details["2D-3"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV overlaps the only (last) exon of {gene_symbol} and No established pathogenic small variants have been report in this exon. \n "

							else:
								if cnv.end >= exons[0][1] and cnv.end <= exons[1][0]:
									if snvs_str:
										# snv_count = 0
										# for snv in snv_list:
										# 	if int(snv[1]) <= cnv.end:
										# 		snv_count += 1
										# ## 2D-2
										# if snv_count:
										rules[f'2D-2_HI-{HI_score}'] = True
										details["2D-2"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV overlaps the last exon of {gene_symbol} only and {snvs_str} established pathogenic small variants have been report in this exon incluing {snvs_str}. \n "

									else:
										## 2D-3
										rules[f'2D-3_HI-{HI_score}'] = True
										details["2D-3"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV overlaps the last exon of {gene_symbol} only and No established pathogenic small variants have been report in this exon. \n "

								## 2D-4
								if cnv.start <= exons[1][1]:
									rules[f'2D-4_HI-{HI_score}'] = True				
									details["2D-4"] += f"* {gene_symbol}: \n ClinGen Haploinsufficiency score: {HI_score}. \n Strand location: minus (-) strand. \n This CNV spans the 3' UTR and includes other exons in addition to the last exon. Nonsense-mediated decay expected to occur. \n "

			elif HI_score == 40:
				if overlap_mode in ["<", "="]:
					rules['2F'] = True
					if "* Genes: \n" not in details["2F"]:
						details["2F"] += "* Genes: \n"
					details["2F"] += f" - {gene_symbol}: ClinGen Haploinsufficiency score: {HI_score}.\n "

			## 2E 基因内缺失
			if HI_score in [3, 2] and overlap_mode == "<":
				tx_refseq_name = row["RefSeq_tx_id"].split("|")[0]
				pvs1_strength, pvs1_class = PVS1_Strength(cnv, tx_refseq_name)
				rules[f'2E-{pvs1_class}_HI-{HI_score}'] = True
				details[f'2E-{pvs1_class}_HI-{HI_score}'] = f""" PVS1 strength is {pvs1_strength} according to ClinGen SVI working group PVS1 specifications.
* ClinGen Haploinsufficiency score: {HI_score}.
* The MANE Select transcript ID: {tx_refseq_name}.
* gnomAD PLI score: {pLI_score}.
* the upper bound of the observed/expected confidence interval: {oe_lof_upper}).
* DECIPHER HI score: {HI_index}.
"""
		#### DEL ####



	## 评估 ClinGen 的致病性/良性区域
	for region in ClinGen_region(cnv, DB_ClinGen_regions):

		region_lst, score, overlap, coverage, mode, time, id = region
		region_chr, region_st, region_fl = region_lst
		region_coord = f"{region_chr}:{region_st}-{region_fl}:{score}({mode})({id})({time})"
		out_record.ClinGen_regions.append(region_coord)

		if "DUP" in cnv.type:
			## 2A-B 与确定为TS的ClinGen region重叠
			if score in [3, 2]:
				#if score == 2 and int(time.split("-")[0]) < 2019:
				#	continue

				if mode in [">", "="]:
					## 2A: 已知TS区域完全重叠
					rules[f'2A_TS-{score}'] = True
					if "* Regions: \n" not in details["2A"]:
						details["2A"] += "* Regions: \n"
					details["2A"] += f" - {id}: ClinGen Triplosensitivity score: {score}. The last evaluated date is {time}\n "

				else:
					region = CNVinfo(region_chr, region_st, region_fl)
					region_genes = count_overlap_genes(region, DB_gene_info)
					cnv_genes = count_overlap_genes(cnv, DB_gene_info)
					
					for gene, gene_info in cnv_genes.items():
						_, cnv_gene_coverage, _, TS_score, _ = gene_info
						if gene not in region_genes or cnv_gene_coverage != 1:
							continue
						if TS_score in ["3", "2"]:
							rules[f'2A_TS-{score}'] = True
							break
						else:
							item = rules.get(f'2A_TS-{score}', 0)
							if type(item) is bool:
								continue
							if coverage < item:
								continue
							rules[f'2A_TS-{score}'] = coverage
											
					if "* Regions: \n" not in details["2A"]:
						details["2A"] += "* Regions: \n"
					details["2A"] += f" - {id}: ClinGen Triplosensitivity score: {score}. The last evaluated date is {time}\n "
					
					
			## 2C-2G 与良性的ClinGen region重叠
			elif score == 40:
				if mode == "=":
					rules['2C'] = True

					if "* Regions: \n" not in details["2C"]:
						details["2C"] += "* Regions: \n"
					details["2C"] += f" - {id}: ClinGen Triplosensitivity score: {score}. The last evaluated date is {time}\n "
				
				elif mode == "<":
					if not cnv_gene_partial_overlaps_count:
						rules['2D'] = True

						if "* Regions: \n" not in details["2D"]:
							details["2D"] += "* Regions: \n"
						details["2D"] += f" - {id}: ClinGen Triplosensitivity score: {score}. The last evaluated date is {time}\n "

					else:
						rules['2E'] = True

						if "* Regions: \n" not in details["2E"]:
							details["2E"] += "* Regions: \n"				
						details["2E"] += f" - {id}: ClinGen Triplosensitivity score: {score}. The last evaluated date is {time}\n "

				else:
					region = CNVinfo(region_chr, region_st, region_fl)
					region_genes = count_overlap_genes(region, DB_gene_info)
					cnv_genes = count_overlap_genes(cnv, DB_gene_info) 
									   
					flag = 0
					for gene, gene_info in cnv_genes.items():
						_, cnv_gene_coverage, _, TS_score, _ = gene_info
						if gene not in region_genes:
							flag = 1
							break
						_, region_gene_coverage, *_ = region_genes[gene] 
						if cnv_gene_coverage > region_gene_coverage and TS_score != "40":
							flag = 1
							break
					
					if not flag:
						rules['2F'] = True
						if "* Regions: \n" not in details["2F"]:
							details["2F"] += "* Regions: \n"				
						details["2F"] += f" - {id}: ClinGen Triplosensitivity score: {score}. The last evaluated date is {time}\n "

					else:
						rules['2G'] = True
						if "* Regions: \n" not in details["2G"]:
							details["2G"] += "* Regions: \n"				
						details["2G"] += f" - {id}: ClinGen Triplosensitivity score: {score}. The last evaluated date is {time}\n "

		#### 如果是 DEL ####
		else:
			## 2A-B 与确定为HI的ClinGen region重叠
			if score in [3, 2]:
				#if score == 2 and int(time.split("-")[0]) < 2019:
				#    continue

				if mode in [">", "="]:
					## 2A: 已知HI区域完全重叠
					rules[f'2A_HI-{score}'] = True

					if "* Regions: \n" not in details["2A"]:
						details["2A"] += "* Regions: \n"				
					details["2A"] += f" - {id}: ClinGen Haploinsufficiency score: {score}. The last evaluated date is {time}\n "
					
				else:
					region = CNVinfo(region_chr, region_st, region_fl)
					region_genes = count_overlap_genes(region, DB_gene_info)
					cnv_genes = count_overlap_genes(cnv, DB_gene_info)

					for gene, gene_info in cnv_genes.items():
						_, gene_coverage, HI_score, _, gene_coord = gene_info
						if gene not in region_genes or gene_coverage != 1:
							continue

						if HI_score not in ["3", "2"]:
							null_variant_nums = gene_null_variant(gene_coord, DB_Clinvar_null_variant)
							# if null_variant_nums >= p_LOF_threshold:
							# 	HI_score = "3"
							if null_variant_nums >= 9:
								HI_score = "3"
							elif null_variant_nums >= 7 and null_variant_nums < 9:
								HI_score = "2"
							elif null_variant_nums >= 5 and null_variant_nums < 7:
								HI_score = "1"
							else:
								HI_score = "0"
							out_record.p_LOF_genes.append(f"{gene}(null variants: {null_variant_nums})")

						if HI_score in ["3", "2"]:
							rules[f'2A_HI-{score}'] = True
							break
						else:
							item = rules.get(f'2A_HI-{score}', 0)
							if type(item) is bool:
								continue
							if coverage < item:
								continue
							rules[f'2A_HI-{score}'] = coverage
						
					if "* Regions: \n" not in details["2A"]:
						details["2A"] += "* Regions: \n"				
					details["2A"] += f" - {id}: ClinGen Haploinsufficiency score: {score}. The last evaluated date is {time}\n "

			## 2F 与良性的ClinGen region重叠
			elif score == 40:
				if mode in ["=", "<"]:
					rules['2F'] = True
					if "* Regions: \n" not in details["2F"]:
						details["2F"] += "* Regions: \n"				
					details["2F"] += f" - {id}: ClinGen Haploinsufficiency score: {score}. The last evaluated date is {time}\n "
					
				else:
					region = CNVinfo(region_chr, region_st, region_fl)
					region_genes = count_overlap_genes(region, DB_gene_info)
					cnv_genes = count_overlap_genes(cnv, DB_gene_info)
									   
					flag = 0
					for gene, gene_info in cnv_genes.items():
						_, cnv_gene_coverage, HI_score, *_ = gene_info
						if gene not in region_genes:
							flag = 1
							break                            
						_, region_gene_coverage, *_ = region_genes[gene]  
						if cnv_gene_coverage > region_gene_coverage and HI_score != "40":
							flag = 1
							break                    

					if not flag:
						rules['2F'] = True

						if "* Regions: \n" not in details["2F"]:
							details["2F"] += "* Regions: \n"				
						details["2F"] += f" - {id}: ClinGen Haploinsufficiency score: {score}. The last evaluated date is {time}\n "
					
					else:
						rules['2G'] = True

						if "* Regions: \n" not in details["2G"]:
							details["2G"] += "* Regions: \n"				
						details["2G"] += f" - {id}: ClinGen Haploinsufficiency score: {score}. The last evaluated date is {time}\n "
	
	## 第3部分评分
	nums = count_procoding_genes(gene_list)
	if "DUP" in cnv.type:
		if nums < 35:
			rules['3A'] = True
			details["3A"] = f"* Number of overlapping protein-coding genes (not from the same family) is {nums}.\n "
		elif nums >= 35 and nums <= 49:
			rules['3B'] = True
			details["3B"] = f"* Number of overlapping protein-coding genes (not from the same family) is {nums}.\n "
		else:
			rules['3C'] = True
			details["3C"] = f"* Number of overlapping protein-coding genes (not from the same family) is {nums}.\n "
	else:
		if nums < 25:
			rules['3A'] = True
			details["3A"] = f"* Number of overlapping protein-coding genes (not from the same family) is {nums}.\n "
		elif nums >= 25 and nums <= 34:
			rules['3B'] = True
			details["3B"] = f"* Number of overlapping protein-coding genes (not from the same family) is {nums}.\n "
		else:
			rules['3C'] = True
			details["3C"] = f"* Number of overlapping protein-coding genes (not from the same family) is {nums}.\n "


	### 4L
	#score = 0
	#dbvar_pathogenoc_set = count_dbvar_pathogenoc_CNV(cnv, DB_dbvar_pathogenic)
	#if dbvar_pathogenoc_set:
	#    details["4L"] += f" Varinat ID:\n "
	#for key, value in dbvar_pathogenoc_set.items():
	#    overlap, coverage, id = value
	#    if overlap >= 0.8 and coverage >= 0.8:
	#        score += 0.05
	#    elif overlap >= 0.5 and coverage >= 0.5:
	#        score += 0.025
	#    else:
	#        continue
	#    out_record.dbvar_pathogenic.append(key)
	#    details["4L"] += f" - {id}\n "

	#score = round(score, 2)
	#if score > 0.15:
	#    score = 0.15
	#if score > 0:
	#    rules['4L'] = score

	### 4N
	#score = 0
	#dbvar_common_set = count_dbvar_common_CNV(cnv, DB_dbvar_common)
	#if dbvar_common_set:
	#    details["4N"] += f" Varinat ID:\n "
	#for key, value in dbvar_common_set.items():
	#    overlap, coverage, var_study = value
	#    if overlap >= 0.8 and coverage >= 0.8:
	#        score += -0.05
	#    elif overlap >= 0.5 and coverage >= 0.5:
	#        score += -0.025
	#    else:
	#        continue
	#    out_record.dbvar_common.append(key)
	#    details["4N"] += f" - {var_study}\n "

	#score = round(score, 2)
	#if score < -0.15:
	#    score = -0.15
	#if score < 0:
	#    rules['4N'] = score
		

	## 4L/LN
	score_p = 0
	score_b = 0
	clinvar_cnv_set = count_ClinVar_CNV(cnv, DB_ClinVar_CNV)
	if clinvar_cnv_set:
		details["4L"] += f" Varinat ID:\n "
		details["4N"] += f" Varinat ID:\n "
	for key, value in clinvar_cnv_set.items():
		overlap, coverage, var_class, var_id = value  
		if any(my_class in var_class for my_class in ["Likely pathogenic", "Pathogenic"]):
			out_record.ClinVar_pathogenic.append(key)
			details["4L"] += f" - {var_id}\n "
			if overlap >= 0.8 and coverage >= 0.8:
				score_p += 0.05
			else:
				score_p += 0.025
		else:
			out_record.ClinVar_benign.append(key)
			details["4N"] += f" - {var_id}\n "
			if overlap >= 0.8 and coverage >= 0.8:
				score_b += -0.05
			else:
				score_b += -0.025
		
	score_p = round(score_p, 2)
	score_b = round(score_b, 2)
	if score_p > 0:
		rules['4L'] = min(score_p, 0.15)
	
	if score_b < 0:
		rules['4N'] = max(score_b, -0.15)

	  
	## 4O
	benign_CNV_set = benign_CNV(cnv, DB_DGV_GS_CNV, DB_DDD_CNV, DB_gnomAD_SV)
	score0 = 0
	score = 0
	cnv_gene_dict = count_overlap_genes(cnv, DB_gene_info)
	for key, cnv_set in benign_CNV_set.items():
		if not cnv_set:
			continue
		for common_cnv in cnv_set:
			coord, mode, overlap, coverage = common_cnv            
			if mode in ["=", "<"]:
				score = -1
			else:
				if overlap < 0.5 or coverage < 0.5:
					continue
				coord_lst = coord.split("(")[0].replace(":", "-").split("-")
				region = CNVinfo(coord_lst[0], coord_lst[1], coord_lst[2])
				region_gene_dict = count_overlap_genes(region, DB_gene_info)
				if any(gene not in region_gene_dict for gene in cnv_gene_dict):
					continue
				
				score = -1
				for gene in cnv_gene_dict:
					_, cnv_coverage, *_ = cnv_gene_dict[gene]
					_, region_coverage, *_ = region_gene_dict[gene]
					#if region_coverage == 1 and cnv_coverage == 1:
					if region_coverage == 1:
						continue
					if cnv_coverage > region_coverage:
						score = score + (cnv_coverage - region_coverage) / 10
			
			if score < score0:
				score0 = score

			overlap = f"{overlap:.2g}"
			coverage = f"{coverage:.2g}"
			if key == "DGV":
				if "* DGV: \n" not in details["4O"]:
					details["4O"] += "* DGV: \n"				
				details["4O"] += f" - {coord}. overlaped percentage of CNV: {overlap}. overlaped percentage of common variant: {coverage}\n "
				out_record.DGV.append(coord)
			
			elif key == "DDD":
				if "* DDD: \n" not in details["4O"]:
					details["4O"] += "* DDD: \n"				
				details["4O"] += f" - {coord}. overlaped percentage of CNV: {overlap}. overlaped percentage of common variant: {coverage}\n "				
				out_record.DECIPHER.append(coord)
			
			else:
				if "* gnomAD: \n" not in details["4O"]:
					details["4O"] += "* gnomAD: \n"				
				details["4O"] += f" - {coord}. overlaped percentage of CNV: {overlap}. overlaped percentage of common variant: {coverage}\n "				
				out_record.gnomAD.append(coord)

	if score0:
		rules['4O'] = score0


	## 评估致病性
	total_points, pathogenicity, score_table= classify(cnv.type, rules)
	out_record.pathogenicity = pathogenicity
	out_record.total_points = total_points
	return score_table, details, out_record


def out_file_helper(CNV, separator):
	score_table, _, out_record = annotation(CNV)
	out_lst = CNV[:-1] + out_record.record_list()
	for value in score_table.values():
		out_lst.append(str(value))
	annotation_line = f"{separator}".join(out_lst)

	return annotation_line


def interpret(input, build, out_file):
	out_file_format = out_file.split(".")[-1]
	if out_file_format in ["tsv", "txt"]:
		separator = "\t"
	elif out_file_format in ["csv"]:
		separator = ","
	else:
		print("Output file format not supported. Please use tsv, txt, or csv.")
		exit(1)

	start_time = time.perf_counter()	
	out_dir = os.path.dirname(out_file)
	
	# 检测输入文件格式并转换
	if any(input.endswith(suffix) for suffix in [".vcf", ".bed"]):		
		if (input).endswith(".vcf"):
			tmp_bed_file_name = os.path.basename(input).replace(".vcf", ".bed")
			tmp_bed_file_path = os.path.join(out_dir, tmp_bed_file_name)
			vcf2bed(input, tmp_bed_file_path)
			input = tmp_bed_file_path
		elif (input).endswith(".bed"):
			input = input
		else:
			print("please input right format!")
			exit(1)
	else:
		CNV_record_file = os.path.join(out_dir, f"{input}.bed")
		with open(CNV_record_file, "w") as f:
			cnv_list = CNV_record_normalize(input)
			print("\t".join(["#chr", "start", "end", "type"]), file = f)
			print("\t".join(cnv_list), file = f)
		input = CNV_record_file

	# 参考基因组版本转换
	if build == "T2T":
		build = "hg38"
		out_dir = os.path.dirname(out_file)
		input = file_build_convert(input, out_dir, chain_file_T2T)

	# 创建CNV对象
	with open(input, "r") as f:
		total_lines = sum(1 for l in f if not l.startswith("#"))

	with open(input, "r") as f, open(out_file, "w") as out:
		raw_header_list = next(f).strip().split("\t")
		header_line = make_header(raw_header_list, separator)
		print(header_line, file = out)
		for line in tqdm(f, total = total_lines, desc='Processing CNVs', unit='CNV'):
			if line.startswith("#"):
				continue
			cnv = line.rstrip().split("\t")
			cnv.append(build)
			result = out_file_helper(cnv, separator)
			print(result, file = out)
	
	end_time = time.perf_counter()	
	print('Elapsed time of annotation and interpretation:', '{0:.1f}'.format(end_time - start_time), 'seconds')
