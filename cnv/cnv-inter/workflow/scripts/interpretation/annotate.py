from interpretation.utils import check_overlap


## 确定真实的蛋白编码基因数量
def count_genes(cnv, db):
	result = db.overlap(cnv.chr, cnv.start, cnv.end)
	nums = 0
	gene_list = []
	for record, overlap, coverage in result:
		record = record._asdict()
		gene_name = record["gene_symbol"].upper()
		gene_list.append(gene_name)
		nums += 1
	gene_list.sort()
	gene_str = ", ".join(gene_list)
	return nums, gene_str



def count_overlap_genes(cnv, db):
	result = db.overlap(cnv.chr, cnv.start, cnv.end)
	gene_dict = {}
	for record, overlap, coverage in result:
		row = record._asdict()
		gene_chr = row["gene_chr"]
		gene_st = int(row["gene_start"]) + 1
		gene_fl = row["gene_end"]
		gene_symbol = row["gene_symbol"]
		HI_score = row["HI_Score"]
		TS_score = row["TS_Score"]
		gene_dict[gene_symbol] = [overlap, coverage, HI_score, TS_score, [gene_chr, gene_st, gene_fl]]
	return gene_dict



## CNV 区带信息
def CNV_cytoband(cnv, db):
	cytobands = db.overlap(cnv.chr, cnv.start, cnv.end)
	cytoband_lst = []
	for record, overlap, coverage in cytobands:
		record = record._asdict()
		chr = record["chr"]
		if chr.startswith("chr"):
			chr = chr.split("chr")[-1]
		cytoband = chr + record["cytoband"]
		cytoband_lst.append(cytoband)
	
	if not cytoband_lst:
		print(cnv.chr, cnv.start, cnv.end)
		return
		
	if len(cytoband_lst) == 1:
		return f"{cytoband_lst[0]}({cytoband_lst[0]})"
	else:
		return f"{cytoband_lst[0]}-{cytoband_lst[-1]}({';'.join(cytoband_lst)})"



## 收集ClinGen中收录的region
def ClinGen_region(cnv, db):
	col_index = "HI_Score" if cnv.type == "DEL" else "TS_Score"
	result = db.overlap(cnv.chr, cnv.start, cnv.end)
	for record, overlap, coverage in result:
		record = record._asdict()
		region_chr = record["chr"]
		region_st = record["start"]
		region_fl = record["end"]
		id = record["ISCA_ID"]
		time = record["date_last_evaluated"]
		clingen_region = [region_chr, region_st, region_fl]
		mode = check_overlap(overlap, coverage) 
		score = record[col_index]
		
		if score in ["Not yet evaluated", "-"]:
			continue

		score = int(score)
		if score not in [3, 2, 40]:
			continue
		
		region = [clingen_region, score, overlap, coverage, mode, time, id]
		yield region



## 获得某个基因的详细信息
def gene_content(result, db):

	gene_info = []

	chr = result["gene_chr"]
	strand = result["strand"]

	exon_st = map(lambda x: int(x) + 1, result["exon_starts"].strip(',').split(','))
	exon_fl = map(int, result["exon_ends"].strip(',').split(','))
	exons = list(zip(exon_st, exon_fl))

	cds = [int(result["CDS_starts"])+1, int(result["CDS_ends"])]
	
	if strand == "+":
		five_prime_UTR =  [int(result["transcript_start"])+1, cds[0]]
		three_prime_UTR = [int(result["transcript_end"])+1, cds[1]]
	else:
		five_prime_UTR =  [int(result["transcript_end"])+1, cds[1]]
		three_prime_UTR = [int(result["transcript_start"])+1, cds[0]]	
	
	gene_chr = result["gene_chr"]
	gene_st = result["gene_start"]
	gene_fl = result["gene_end"]
	result_variants = db.overlap(gene_chr, gene_st, gene_fl)
	snv_list = []
	for record, overlap, coverage in result_variants:
		record = record._asdict()
		snv_str = "-".join([record["chr"], str(record["start"]), str(record["end"]), f'{record["HGVS"]}-{record["stat"]}-{record["classification"]}'])
		snv_list.append(snv_str)
	snvs_str = ", ".join(snv_list)

	gene_info.append(strand)
	gene_info.append(five_prime_UTR)
	gene_info.append(three_prime_UTR)
	gene_info.append(cds)
	gene_info.append(exons)
	gene_info.append(snvs_str)

	return gene_info


## 获得 null variant 数量
def gene_null_variant(gene_coord, db):
	result = db.overlap(gene_coord[0], gene_coord[1], gene_coord[2])
	nums = len(list(result))
	return nums


def benign_CNV(cnv, DB_DGV, DB_DDD, DB_gnomAD):

	common_cnv_list = {
		"DGV": [],
		"DDD": [],
		"gnomAD": []
	}

	## DGV
	dgv = DB_DGV.overlap(cnv.chr, cnv.start, cnv.end)
	for record, overlap, coverage in dgv:
		record = record._asdict()
		DGV_region = [record["chr"], record["start"], record["end"]]
		DGV_type = record["type"]
		if DGV_type != cnv.type:
			continue
		id = record["ID"]
		af = record["frequency"]
		if eval(af) < 0.01:
			continue
		af = f"{eval(af):.2g}"
		overlap_mode = check_overlap(overlap, coverage)
		coord = f"{DGV_region[0]}:{DGV_region[1]}-{DGV_region[2]}-{id}(AF={af})"
		if [coord, overlap_mode, overlap, coverage] not in common_cnv_list["DGV"]:
			common_cnv_list["DGV"].append([coord, overlap_mode, overlap, coverage])

	## DDD
	col_index = "del_freq" if cnv.type == "DEL" else "dup_freq"
	decipher = DB_DDD.overlap(cnv.chr, cnv.start, cnv.end)
	for record, overlap, coverage in decipher:
		record = record._asdict()
		DDD_region = [record["chr"], record["start"], record["end"]]
		af = record[col_index]
		if eval(af) < 0.01:
			continue
		af = f"{eval(af):.2g}"
		overlap_mode = check_overlap(overlap, coverage)
		coord = f"{DDD_region[0]}:{DDD_region[1]}-{DDD_region[2]}(AF={af})"
		if [coord, overlap_mode, overlap, coverage] not in common_cnv_list["DDD"]:
			common_cnv_list["DDD"].append([coord, overlap_mode, overlap, coverage])
	
	## gnomAD	
	gnomad = DB_gnomAD.overlap(cnv.chr, cnv.start, cnv.end)
	for record, overlap, coverage in gnomad:
		record = record._asdict()
		gnomAD_region = [record["chr"], record["start"], record["end"]]
  
		gnomad_type = record["svtype"]
		if gnomad_type != cnv.type:
			continue

		if cnv.build == "hg38":
			af_columns = [col for col in record.keys() if col.startswith('AF')]
		else:
			af_columns = [col for col in record.keys() if col.endswith('AF')]	

		if all(eval(record[AF]) < 0.01 for AF in af_columns):
			continue
		af_list = [f"{AF}={eval(record[AF]): .2g}" for AF in af_columns]
		af_str = ";".join(af_list)

		overlap_mode = check_overlap(overlap, coverage)
		coord = f"{gnomAD_region[0]}:{gnomAD_region[1]}-{gnomAD_region[2]}({af_str})"
		if [coord, overlap_mode, overlap, coverage] not in common_cnv_list["gnomAD"]:
			common_cnv_list["gnomAD"].append([coord, overlap_mode, overlap, coverage])		

	return common_cnv_list


def count_dbvar_pathogenoc_CNV(cnv, db):
	result = db.overlap(cnv.chr, cnv.start, cnv.end)
	CNV_lst = {}
	#type_dict = {
	#    "DEL" : ["copy_number_loss", "deletion", "copy_number_loss;deletion"],
	#    "DUP" : ["copy_number_gain", "duplication", "copy_number_gain;duplication"],
	#    "TDUP" : ["copy_number_gain", "duplication", "copy_number_gain;duplication"]
	#}
	for record, overlap, coverage in result:
		record = record._asdict()
		chr = record["chr"]
		start = record["outermost_start"]
		end = record["outermost_stop"]
		variant_type = record["variant_type"]
		clinical_assertion = record["clinical_assertion"]
		id = record["variant"]
		if variant_type not in cnv.type:
			continue
		CNV_record = f"{chr}:{start}-{end}({id})"
		CNV_lst[CNV_record] = [overlap, coverage, id]
	return CNV_lst


def count_dbvar_common_CNV(cnv, db):
	result = db.overlap(cnv.chr, cnv.start, cnv.end)
	CNV_lst = {}
	#type_dict = {
	#    "DEL" : ["copy_number_loss", "deletion", "copy_number_loss;deletion"],
	#    "DUP" : ["copy_number_gain", "duplication", "copy_number_gain;duplication"],
	#    "TDUP" : ["copy_number_gain", "duplication", "copy_number_gain;duplication"]
	#}
	for record, overlap, coverage in result:
		record = record._asdict()
		chr = record["chr"]
		start = record["outermost_start"]
		end = record["outermost_stop"]
		variant_type = record["variant_type"]
		var_study = record["study"]
		id = record["variant"]
		if variant_type not in cnv.type:
			continue
		CNV_record = f"{chr}:{start}-{end}({var_study})"
		CNV_lst[CNV_record] = [overlap, coverage, var_study]
	return CNV_lst


def count_ClinVar_CNV(cnv, db):
	result = db.overlap(cnv.chr, cnv.start, cnv.end)
	CNV_lst = {}
	for record, overlap, coverage in result:
		if overlap < 0.5 or coverage < 0.5:
			continue
		record = record._asdict()
		chr = record["chr"]
		start = record["start"]
		end = record["end"]
		cnv_type = record["type"]
		var_class = record["variant_class"]
		var_id = record["variant_ID"]
		var_accession_ID = record["accession_ID"]
		#status = record["status"]
		if cnv_type not in cnv.type:
			continue
		CNV_record = f"{chr}:{start}-{end}({var_id};{var_accession_ID})"
		CNV_lst[CNV_record] = [overlap, coverage, var_class, var_id]
	return CNV_lst
