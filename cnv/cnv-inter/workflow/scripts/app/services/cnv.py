import os
import traceback
from typing import List
import requests

from loguru import logger
from app.schemas.msg import ServiceError

from interpretation.interpret import annotation
from interpretation.interpret import interpret
from interpretation.utils import CNV_record_normalize
from interpretation.utils import CNV_record_build_convert
from interpretation.resources import chain_file_T2T


# from app.schemas.acmg import EvidenceParameter
# from app.schemas.msg import ServiceError

		
def interpret_cnv_record(cnv_str: str, genome_build: str = "hg38"):
	"""
	cnv_str example: "chr1:10000-20000-DEL"
	"""

	res_dict = {"total_score": "",
		"pathogenicity": "",
		"retrieval_records": {},
		"base_info": {
			"type": "",
			"coord": "",
			"cytoband": "",
			"size": "",
			"interitance": ""
		},
		"rules_table": {}
	}

	# FILE_DIR = "/nas/CNVSeeker/cnvseeker/cnv/cnv-inter/workflow/scripts/"
	FILE_DIR = os.getcwd()
	tmp_dir = FILE_DIR + "/web_output"

	cnv_record_list = CNV_record_normalize(cnv_str)
	chr, st, fl, cnv_type = cnv_record_list

	if genome_build == "T2T":
		cnv_record_list = CNV_record_build_convert(cnv_record_list, tmp_dir, chain_file_T2T)
	cnv_record_list.append(genome_build)
	score_table, details, out_record = annotation(cnv_record_list)

	if cnv_type == "DEL":
		sub_type = "Copy number loss"
	else:
		sub_type = "Copy number gain"

	coord = f"{chr}:{st}-{fl}-{cnv_type} ({genome_build})"
	res_dict["total_score"] = out_record.total_points
	res_dict["pathogenicity"] = out_record.pathogenicity
	res_dict["base_info"]["type"] = sub_type
	res_dict["base_info"]["coord"] = coord

	cytoband = out_record.cytoband.split("(")[0]
	res_dict["base_info"]["cytoband"] = cytoband
	cytoband_lst = out_record.cytoband.split("(")[1].split(")")[0].split(";")

	keywords_dict = {
		"logic": "OR",
		"conditions": [
			{"logic": "AND", "conditions": [{"logic": "OR", "words": []}, {"logic": "OR", "words": []}]},
			{"logic": "OR", "words": []}
		]
	}
	type_list, cytoband_words_lst_1 = retrieval_records_1(cytoband_lst, cnv_type)
	cytoband_words_lst_2 = retrieval_records_2(cytoband_lst, chr, cnv_type)
	keywords_dict["conditions"][0]["conditions"][0]["words"] = type_list
	keywords_dict["conditions"][0]["conditions"][1]["words"] = cytoband_words_lst_1
	keywords_dict["conditions"][1]["words"] = cytoband_words_lst_2
	res_dict["retrieval_records"] = keywords_dict

	res_dict["base_info"]["size"] = out_record.size

	for rule, score_record in score_table.items():
				
		if score_record:
			for record in score_record.split(";"):
				if "_" in score_record:
					rule = score_record.split("_")[0]
				else:
					rule = score_record.split("=")[0]
				score = score_record.split("=")[-1]
				detail = details.get(rule, None)
				show = True
				if rule not in res_dict["rules_table"].keys():
					res_dict["rules_table"][rule] = {"score" : score, "detail": detail, "show": show}
					continue
				if (eval(score) > 0 and eval(score) > eval(res_dict["rules_table"][rule]["score"])) or (eval(score) < 0 and eval(score) < eval(res_dict["rules_table"][rule]["score"])):
					res_dict["rules_table"][rule]["score"] = score
		else:
			score = ""		
			detail = details.get(rule, None)
			show = False
			res_dict["rules_table"][rule] = {"score" : score, "detail": detail, "show": show}
	
	return res_dict


def interpret_pedigree_cnv_record(cnv_str: str, genome_build: str = "hg38", relation_list: list = [{"relation": "father", "path": None}, {"relation": "mother", "path": None}]):
	"""
	cnv_str example: "chr1:10000-20000-DEL"
	"""

	res_dict = {"total_score": "",
		"pathogenicity": "",
		"retrieval_records": {},
		"base_info": {
			"type": "",
			"coord": "",
			"cytoband": "",
			"size": ""
		},
		"rules_table": {}
	}

	# FILE_DIR = "/nas/CNVSeeker/cnvseeker/cnv/cnv-inter/workflow/scripts/"
	FILE_DIR = os.getcwd()
	tmp_dir = FILE_DIR + "/web_output"

	cnv_record_list = CNV_record_normalize(cnv_str)
	chr, st, fl, cnv_type = cnv_record_list

	if genome_build == "T2T":
		cnv_record_list = CNV_record_build_convert(cnv_record_list, tmp_dir, chain_file_T2T)
	cnv_record_list.append(genome_build)
	score_table, details, out_record = annotation(cnv_record_list)

	if cnv_type == "DEL":
		sub_type = "Copy number loss"
	else:
		sub_type = "Copy number gain"

	coord = f"{chr}:{st}-{fl}-{cnv_type}"
	res_dict["total_score"] = out_record.total_points
	res_dict["pathogenicity"] = out_record.pathogenicity
	res_dict["base_info"]["type"] = sub_type
	res_dict["base_info"]["coord"] = coord

	cytoband = out_record.cytoband.split("(")[0]
	res_dict["base_info"]["cytoband"] = cytoband
	cytoband_lst = out_record.cytoband.split("(")[1].split(")")[0].split(";")

	keywords_dict = {
		"logic": "OR",
		"conditions": [
			{"logic": "AND", "conditions": [{"logic": "OR", "words": []}, {"logic": "OR", "words": []}]},
			{"logic": "OR", "words": []}
		]
	}
	type_list, cytoband_words_lst_1 = retrieval_records_1(cytoband_lst, cnv_type)
	cytoband_words_lst_2 = retrieval_records_2(cytoband_lst, chr, cnv_type)
	keywords_dict["conditions"][0]["conditions"][0]["words"] = type_list
	keywords_dict["conditions"][0]["conditions"][1]["words"] = cytoband_words_lst_1
	keywords_dict["conditions"][1]["words"] = cytoband_words_lst_2
	res_dict["retrieval_records"] = keywords_dict

	res_dict["base_info"]["size"] = out_record.size

	for rule, score_record in score_table.items():
		if score_record:
			score = eval(score_record.split("=")[-1])
			if (rule == "2I" and cnv_type == "DUP") or (rule in ["2C", "2D", "2E"] and cnv_type == "DEL"):
				rule = score_record.split("=")[0]
			detail = details.get(rule, None)
			show = True
		else:
			score = score_record
			detail = details.get(rule, None)
			show = False		
		res_dict["rules_table"][rule] = {"score" : score, "detail": detail, "show": show}
	
	interitance = is_denovo(cnv_record_list, relation_list)
	res_dict["base_info"]["interitance"] = interitance

	return res_dict


def interpret_cnv_records(cnv_list: List[str], outdir: str, outname: str, genome_build:str='hg38', callback:str=None, is_backgroud: bool = True):
	
	try:
		outdir = outdir if outdir else os.getcwd() + "/web_output"
		outname = outname if outname else "result.tsv"

		if not os.path.isdir(outdir):
			os.makedirs(outdir)
			logger.debug(f"Create directory: <{outdir}>.")
		
		infile = os.path.join(outdir, "cnv_list.bed")
		with open(infile, "w") as f:
			for cnv_str in cnv_list:
				cnv_record_list = CNV_record_normalize(cnv_str)	
				print("\t".join(cnv_record_list), file = f)

		outfile = os.path.join(outdir, outname)

		if os.path.isfile(infile):
			logger.warning(f"The file <{infile}> already exists!")
		else:
			logger.info(f"Output file is <{infile}>.")

		if os.path.isfile(outfile):
			logger.warning(f"The file <{outfile}> already exists!")
		else:
			logger.info(f"Output file is <{outfile}>.")

		cores = 20
		interpret(infile, cores, genome_build, outfile)
	
		url = callback + "?code=1"
		response = requests.get(url)
		if response.status_code != 200:
			logger.warning(f"response from {url}: {response.status_code}")

		return outfile
	
	except Exception as e:
		res = ServiceError(
			msg = str(e),
			traceback = traceback.format_exc()
		)
		logger.exception(f"{e}\n")

		if is_backgroud:
			url = callback + "?code=0"
			requests.get(url=url)
		else:
			return res
		

def retrieval_records_1(cytoband_lst, cnv_type):
	
	if cnv_type == "DEL":
		type_list = ["deletion", "copy number loss"]
	else:
		type_list = ["duplication", "copy number gain"]

	cytoband_words_lst = [i for i in cytoband_lst]
	for i, cytoband1 in enumerate(cytoband_lst):
		for cytoband2 in cytoband_lst[i+1:]:
			cytoband_words_lst.append(f"{cytoband1}-{cytoband2}")

	return type_list, cytoband_words_lst


def retrieval_records_2(cytoband_lst, chr, cnv_type):

	if chr.startswith("chr"):
		chr = chr.split("chr")[-1]

	cytoband_lst_new = []
	for cytoband in cytoband_lst:
		if "p" in cytoband:
			cytoband_lst_new.append("p" + cytoband.split("p")[-1])
		elif "q" in cytoband:
			cytoband_lst_new.append("q" + cytoband.split("q")[-1])
		else:
			continue

	cytoband_words_lst = []
	concat_str_lst = [";", "; ", "-", "", " to "]
	for i, cytoband1 in enumerate(cytoband_lst_new):
		cytoband_words_lst.append(cytoband1)
		for cytoband2 in cytoband_lst_new[i+1:]:
			if cytoband2 not in cytoband_words_lst:
				cytoband_words_lst.append(cytoband2)
				for concat_str in concat_str_lst:
					cytoband_words_lst.append(f"{cytoband1}{concat_str}{cytoband2}")

	final_words_lst = []
	for cytoband in cytoband_words_lst:
		final_words_lst.append(f"{cnv_type.lower()}({chr})({cytoband})")

	return final_words_lst

# print(retrieval_records_1(["7p15.2", "7p15.1", "7p14.3", "7p14.2", "7p14.1"], "DUP"))
# print(retrieval_records_2(["7p15.2", "7p15.1", "7p14.3", "7p14.2", "7p14.1"], "chr7", "DUP"))


def is_denovo(cnv, relation_list):
	relation_result_dict = {}
	cnv_str = "-".join(cnv)
	tmp_dir = os.getcwd() + f"/web_output/tmp_{cnv_str}/"
	if not os.path.exists(tmp_dir):
		os.makedirs(tmp_dir)
	
	cnv_tmp_bed = tmp_dir + "cnv.bed"
	with open(cnv_tmp_bed, "w") as f:
		print("\t".join(cnv), file = f)

	for relation_dict in relation_list:
		relation = relation_dict.get("relation")
		file_path = relation_dict.get("path")
		if file_path.endswith(".vcf"):
			relation_bed_path = tmp_dir + relation + ".bed"
			vcf2bed(file_path, relation_bed_path)
			file_path = relation_bed_path
		out_path = tmp_dir + relation + "_result.bed"
		os.system(f"bedtools intersect -a {cnv_tmp_bed} -b {file_path} -wo -f 0.7 -r > {out_path}")
		if os.path.getsize(out_path):
			relation_result_dict[relation] = 1
		else:
			relation_result_dict[relation] = 0

	if relation_result_dict["father"] + relation_result_dict["mother"] == 0:
		return "This CNV is de novo"
	elif relation_result_dict["father"] == 1:
		return "This CNV is interited from father"
	elif relation_result_dict["mother"] == 1:
		return "This CNV is interited from mother"
	else:
		return "This CNV is interited from parents"
	

def vcf2bed(vcf_path, bed_path):
	with open(bed_path, "w") as out, open(vcf_path, "r") as f:
		for line in f:
			if line.startswith("#"):
				continue
			line = line.rstrip().split("\t")
			sv_type = line[7].split("SVTYPE=")[1].split(";")[0]
			if "SVTYPE" in line[7]:
				sv_type = line[7].split("SVTYPE=")[1].split(";")[0]
			else:
				sv_type = line[4][1:-1]

			chr = line[0]
			if not chr.startswith("chr"):
				chr = "chr" + chr
				
			st = str(int(line[1]) - 1)
			info = line[7]
			if "SVLEN=" in info:
				cnv_len = abs(int(info.split("SVLEN=")[1].split(";")[0]))
				fl = str(int(st) + cnv_len)
			elif "END=" in info:
				if info.startswith("END="):
					fl = info.split("END=")[1].split(";")[0]
				else:
					fl = info.split(";END=")[1].split(";")[0]
			else:
				continue

			print("\t".join([chr, st, fl, sv_type]), file = f)
