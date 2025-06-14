from interpretation.autopvs1.cnv import CNVRecord, PVS1CNV
from interpretation.autopvs1.utils import get_transcript
from interpretation.autopvs1.read_data import transcripts_hg38, transcripts_hg19

from multiprocessing import Pool
from tqdm import tqdm
from collections import namedtuple
import pysam
import os
import copy



SUGGESTED_SCORES = {
	'DEL': {
		'1A': 0, '1B': -0.6,
		# HI=3
		'2A_HI-3': 1, 
		'2C-1_HI-3': 0.9, '2C-2_HI-3': 0,
		'2D-1_HI-3': 0, '2D-2_HI-3': 0.9, '2D-3_HI-3': 0.3, '2D-4_HI-3': 0.9,
		'2E-1_HI-3': 0.9, '2E-2_HI-3': 0.45, '2E-3_HI-3': 0.3, '2E-4_HI-3': 0.15, '2E-5_HI-3': 0,   
		# HI=2
		'2A_HI-2': 1, 
		'2C-1_HI-2': 0.9, '2C-2_HI-2': 0,
		'2D-1_HI-2': 0, '2D-2_HI-2': 0.9, '2D-3_HI-2': 0.3, '2D-4_HI-2': 0.9,
		'2E-1_HI-2': 0.9, '2E-2_HI-2': 0.45, '2E-3_HI-2': 0.3, '2E-4_HI-2': 0.15, '2E-5_HI-2': 0,
		# # HI=2
		# '2A_HI-2': 0.75, 
		# '2C-1_HI-2': 0.75, '2C-2_HI-2': 0,
		# '2D-1_HI-2': 0, '2D-2_HI-2': 0.75, '2D-3_HI-2': 0.15, '2D-4_HI-2': 0.75,
		# '2E-1_HI-2': 0.75, '2E-2_HI-2': 0.3, '2E-3_HI-2': 0.15, '2E-4_HI-2': 0.1, '2E-5_HI-2': 0,
		#
		'2B': 0,
		#
		'2F': -1, '2G': 0, '2H': 0.15, 
		'2I': 0, '2J': 0, '2K': 0.45, '2L': 0, # DUP specific
		'3A': 0, '3B': 0.45, '3C': 0.9,
		'4A': 0, '4B': 0, '4C': 0, '4D': 0,
		'4E': 0, '4F': 0, '4G': 0, '4H': 0, '4I': 0, '4J': 0, '4K': 0, '4L': 0, '4M': 0, '4N': 0,
		'4O': -1,
		'5A': 0, '5B': 0, '5C': 0, '5D': 0, '5E': 0, '5F': 0, '5G': 0, '5H': 0
	},
	'DUP': {
		'1A': 0, '1B': -0.6,
		# TS=3
		'2A_TS-3': 1, 
		'2I-1_HI-3': 0.9, '2I-2_HI-3': 0.45, '2I-3_HI-3': 0,         
		# TS=2
		'2A_TS-2': 1,
		'2I-1_HI-2': 0.9, '2I-2_HI-2': 0.45, '2I-3_HI-2': 0,
		# # TS=2
		# '2A_TS-2': 0.75,
		# '2I-1_HI-2': 0.75, '2I-2_HI-2': 0.3, '2I-3_HI-2': 0,
		#
		'2B': 0,
		'2C': -1, '2D': -1, '2E': 0, '2F': -1, '2G': 0, '2H': 0, 
		#
		'2J': 0, '2K': 0.45, '2L': 0,
		'3A': 0, '3B': 0.45, '3C': 0.9,
		'4A': 0, '4B': 0, '4C': 0, '4D': 0, 
		'4E': 0, '4F': 0, '4G': 0, '4H': 0, '4I': 0, '4J': 0, '4K': 0, '4L': 0, '4M': 0, '4N': 0,
		'4O': -1,
		'5A': 0, '5B': 0, '5C': 0, '5D': 0, '5E': 0, '5F': 0, '5G': 0, '5H': 0
	}
}


#SUGGESTED_SCORES = {
#    'DEL': {
#        '1A': 0, '1B': -0.6,
#        # HI=3
#        '2A_HI-3': 1, 
#        '2C-1_HI-3': 0.9, '2C-2_HI-3': 0,
#        '2D-1_HI-3': 0, '2D-2_HI-3': 0.9, '2D-3_HI-3': 0.3, '2D-4_HI-3': 0.9,
#        '2E-1_HI-3': 0.9, '2E-2_HI-3': 0.45, '2E-3_HI-3': 0.3, '2E-4_HI-3': 0.15, '2E-5_HI-3': 0,   
#        # HI=2
#        '2A_HI-2': 0.75, 
#        '2C-1_HI-2': 0.75, '2C-2_HI-2': 0,
#        '2D-1_HI-2': 0, '2D-2_HI-2': 0.75, '2D-3_HI-2': 0.15, '2D-4_HI-2': 0.75,
#        '2E-1_HI-2': 0.75, '2E-2_HI-2': 0.3, '2E-3_HI-2': 0.15, '2E-4_HI-2': 0.1, '2E-5_HI-2': 0,
		
#        #
#        '2B': 0,
#        #
#        '2F': -1, '2G': 0, '2H': 0.15, 
#        '2I': 0, '2J': 0, '2K': 0.45, '2L': 0, # DUP specific
#        '3A': 0, '3B': 0.45, '3C': 0.9,
#        '4A': 0, '4B': 0, '4C': 0, '4D': 0,
#        '4E': 0, '4F': 0, '4G': 0, '4H': 0, '4I': 0, '4J': 0, '4K': 0, '4L': 0, '4M': 0, '4N': 0,
#        '4O': -1,
#        '5A': 0, '5B': 0, '5C': 0, '5D': 0, '5E': 0, '5F': 0, '5G': 0, '5H': 0
#    },
#    'DUP': {
#        '1A': 0, '1B': -0.6,
#        # TS=3
#        '2A_TS-3': 1, 
#        '2I-1_HI-3': 0.9, '2I-2_HI-3': 0.45, '2I-3_HI-3': 0,         
#        # TS=2
#        '2A_TS-2': 0.75,
#        '2I-1_HI-2': 0.75, '2I-2_HI-2': 0.3, '2I-3_HI-2': 0,
#        #
#        '2B': 0,
#        '2C': -1, '2D': -1, '2E': 0, '2F': -1, '2G': 0, '2H': 0, 
#        #
#        '2J': 0, '2K': 0.45, '2L': 0,
#        '3A': 0, '3B': 0.45, '3C': 0.9,
#        '4A': 0, '4B': 0, '4C': 0, '4D': 0, 
#        '4E': 0, '4F': 0, '4G': 0, '4H': 0, '4I': 0, '4J': 0, '4K': 0, '4L': 0, '4M': 0, '4N': 0,
#        '4O': -1,
#        '5A': 0, '5B': 0, '5C': 0, '5D': 0, '5E': 0, '5F': 0, '5G': 0, '5H': 0
#    }
#}


DETAIL_TABLE = {
	"DUP" : {
		"1A" : "",
		"1B" : "",
		"2A" : "Completely overlaped established Triplosensitive genes or regions: \n",
		"2B" : "Partially overlaped established Triplosensitive regions: \n",
		"2C" : "Completely overlaped established benign genes or regions and similar in size: \n",
		"2D" : "Fully contained within the established benign region: \n",
		"2E" : "Fully contained within the established benign region: \n",
		"2F" : "Partially overlaped established benign regions: \n",
		"2G" : "Partially overlaped established benign regions: \n",
		"2H" : "Genes overlapping with the CNV that are predicted to be Haploinsufficient: \n",
		"2I-1" : "",
		"2I-2" : "",
		"2I-3" : "",
		"2J" : "One breakpoint is within an established HI gene: \n",
		"2K" : "One breakpoint is within an established HI gene: \n",
		"2L" : "One breakpoint is within gene(s) of UNKNOWN established clinical significance: \n",
		"4O" : "Overlaped common population variants (MAF > 1%) from DGV, gnomAD, and DECIPHER: \n",
		"4L" : "Overlaped pathogenic/likely pathogenic CNV in dbVar: \n",
		"4N" : "Overlaped benign/likely benign CNV in dbVar: \n"
	},
	"DEL" : {
		"1A" : "",
		"1B" : "",
		"2A" : "Completely overlaped Haploinsufficient genes or regions: \n",
		"2B" : "Partially overlaped established Haploinsufficient regions: \n",
		"2C-1" : "Partial overlaps with the 5' end of an established HI gene: \n",
		"2C-2" : "Partial overlaps with the 5' end of an established HI gene: \n",
		"2D-1" : "Partial overlaps with the 3' end of an established HI gene: \n",
		"2D-2" : "Partial overlaps with the 3' end of an established HI gene: \n",
		"2D-3" : "Partial overlaps with the 3' end of an established HI gene: \n",
		"2D-4" : "Partial overlaps with the 3' end of an established HI gene: \n",	
		"2E-1" : "",
		"2E-2" : "",
		"2E-3" : "",
		"2E-4" : "",
		"2E-5" : "",
		"2F" : "Fully contained within the established benign genes or regions: \n",
		"2G" : "Overlaps established benign regions, but includes additional genomic material: \n",
		"2H" : "Genes overlapping with the CNV that are predicted to be Haploinsufficient: \n",
		"4O" : "Overlaped common population variants (MAF > 1%) from DGV, gnomAD, and DECIPHER: \n",
		"4L" : "Overlaped pathogenic/likely pathogenic CNV in dbVar: \n",
		"4N" : "Overlaped benign/likely benign CNV in dbVar: \n"
	}
}


RULE_TABLE = {
	'1A': 0, '1B': 0,
	'2A': 0, '2B': 0,
	'2C': 0, '2D': 0, '2E': 0, '2F': 0, '2G': 0, '2H': 0, 
	'2I': 0, '2J': 0, '2K': 0, '2L': 0, # DUP specific
	'3A': 0, '3B': 0, '3C': 0,
	'4A': 0, '4B': 0, '4C': 0, '4D': 0, 
	'4E': 0, '4F': 0, '4G': 0, '4H': 0, '4I': 0, '4J': 0, '4K': 0, '4L': 0, '4M': 0, '4N': 0,
	'4O': 0,
	'5A': 0, '5B': 0, '5C': 0, '5D': 0, '5E': 0, '5F': 0, '5G': 0, '5H': 0
}


class DataBase:
	def __init__(self, path):
		self._path = path
		self._tbx = pysam.Tabixfile(path)
		self._fields = namedtuple('Record', self._tbx.header[-1].strip('#').split('\t'))

	def fetch(self, chrom, start, end):
		"""
		查找并生成给定基因组位置的记录
		:param chrom: 染色体编号
		:param start: 起始位置
		:param end: 终止位置
		:return: 记录生成器
		"""
		try:
			for record in self._tbx.fetch(chrom, start, end):
				chrom, start, end, *fields = record.split('\t')
				start, end = int(start), int(end)
				# 所得记录按照表头组装为namedtuple方便使用
				yield self._fields(chrom, start, end, *fields)
		except ValueError:
			yield from ()

	def overlap(self, chrom, start, end):
		"""
		查找并生成给定基因组位置的记录，同时计算两者之间的重叠程度
		:param chrom: 染色体编号
		:param start: 起始位置
		:param end: 终止位置
		:return: 记录与重叠程度生成器
		"""
		length = end - start
		for record in self.fetch(chrom, start, end):
			_, overlap_start, overlap_end, _ = sorted((start, end, record[1], record[2]))
			overlap = overlap_end - overlap_start
			if record[2] == record[1]:
				continue
			yield record, overlap / length, overlap / (record[2] - record[1])


def check_overlap(overlap, coverage):
	if overlap < 1 and coverage < 1:
		overlap_mode = "1"
	elif overlap == 1 and coverage < 1:
		overlap_mode = "<"
	elif overlap < 1 and coverage == 1:
		overlap_mode = ">"
	else:
		overlap_mode = "="

	return overlap_mode


class CNVinfo:
	def __init__(self, chr, start, end, type = None, build = None):
		if "chr" not in chr:
			chr = "chr" + chr
		start = int(start)
		end   = int(end)

		self.chr   = chr
		self.start = start
		self.end   = end
		self.type  = type
		self.build = build

	def size(self):
		length = self.end - self.start + 1
		if length < 1000:
			return f"{float('%.2f' % length)} Bp"
		elif length < 1000000:
			return f"{float('%.2f' % (length / 1000))} Kb"
		else:
			return f"{float('%.2f' % (length / 1000000))} Mb"
	
	def region(self):
		return [self.chr, str(self.start), str(self.end)]


class OutRecord:
	def __init__(self):
		#self.chr             = ""
		#self.start           = 0
		#self.end             = 0
		#self.type            = ""
		self.build           = ""
		self.size			 = ""
		self.cytoband		 = ""
		self.gene_count		 = 0
		self.genes           = list()
		self.HI_genes        = list()
		self.p_HI_genes      = list()
		self.p_LOF_genes	 = list()
		self.pLI_loeuf_HI	 = list()
		self.TS_genes        = list()
		self.OMIM_genes      = list()
		self.ClinGen_regions = list()
		self.DGV             = list()
		self.DECIPHER        = list()
		self.gnomAD          = list()
		self.ClinVar_pathogenic = list()          
		self.ClinVar_benign     = list()
		self.pathogenicity   = "P"
		self.total_points    = "1.0"

	def list_concat(self, lst):
		if lst:
			return ";".join(lst)
		else:
			return "-"

	def list_sort(self, lst):
		lst.sort()
		return lst

	def record_list(self):	
		record_lst = list()
		#record_lst.append(self.chr)
		#record_lst.append(self.start)
		#record_lst.append(self.end)
		#record_lst.append(self.type)
		record_lst.append(self.build)
		record_lst.append(self.size)
		record_lst.append(self.cytoband)
		record_lst.append(self.gene_count)
		record_lst.append(self.list_concat(self.list_sort(self.genes)))
		record_lst.append(self.list_concat(self.list_sort(self.HI_genes)))
		record_lst.append(self.list_concat(self.list_sort(self.p_HI_genes)))
		record_lst.append(self.list_concat(self.list_sort(self.pLI_loeuf_HI)))
		record_lst.append(self.list_concat(self.list_sort(self.p_LOF_genes)))
		record_lst.append(self.list_concat(self.list_sort(self.TS_genes)))
		record_lst.append(self.list_concat(self.list_sort(self.OMIM_genes)))
		record_lst.append(self.list_concat(self.list_sort(self.ClinGen_regions)))
		record_lst.append(self.list_concat(self.list_sort(self.DGV)))
		record_lst.append(self.list_concat(self.list_sort(self.DECIPHER)))
		record_lst.append(self.list_concat(self.list_sort(self.gnomAD)))
		record_lst.append(self.list_concat(self.list_sort(self.ClinVar_pathogenic)))
		record_lst.append(self.list_concat(self.list_sort(self.ClinVar_benign)))
		record_lst.append(self.pathogenicity)
		record_lst.append(self.total_points)
		
		return [str(record) for record in record_lst]


## 对CNV进行分类
def classify(type, rules):
	
	total_points = 0
	pathogenicity = ""
	file_score_table = copy.deepcopy(RULE_TABLE)
	score_table = copy.deepcopy(RULE_TABLE)
			
	## 计算总分
	for rule in rules.keys():
		short_rule = rule[0:2]
		if type == "TDUP":
			type = "DUP"

		rule_score = SUGGESTED_SCORES[type][rule]
		if rules[rule] not in [True, False]:
			rule_score = rules[rule]

		score_record = f"{rule}={rule_score}"

		if (score_table[short_rule] and ((rule_score > 0 and rule_score > score_table[short_rule]) or (rule_score < 0 and rule_score < score_table[short_rule]))) or (not score_table[short_rule]):
			score_table[short_rule] = rule_score

		org_score_record = file_score_table[short_rule]
		if org_score_record == 0:
			new_score_record = score_record
		else:
			new_score_record = ";".join([org_score_record, score_record])
		file_score_table[short_rule] = new_score_record

	total_points = sum(score_table.values())

	## 致病性分类
	if total_points >= 0.99:
		pathogenicity = "Pathogenic (P)"
	elif total_points >= 0.90 and total_points < 0.99:
		pathogenicity = "Likely pathogenic (LP)"
	elif total_points > -0.99 and total_points <= -0.90:
		pathogenicity = "Likely benign (LB)"
	elif total_points <= -0.99:
		pathogenicity = "Benign (B)"
	else:
		pathogenicity = "Variant of uncertain significance (VUS)"
	
	return total_points, pathogenicity, file_score_table


def make_header(raw_header_list, separator):
	#header_list = ["chr", "start", "end", "type", "size (bp)", "CytoBand",
	#            "gene_count", "genes", "known_HI_genes", "predicted_HI_genes", "pLI_loeuf_HI", "known_TS_genes",  "OMIM",
	#            "ClinGen_regions", "DGV", "DECIPHER", "gnomAD", "dbvar_pathogenic", "dbvar_common", 
	#            "pathogenicity", "total_points"]

	header_list = ["build", "size (bp)", "CytoBand",
				"gene_count", "genes", "known_HI_genes", "predicted_HI_genes", "pLI_loeuf_HI", "predicted_LOF_genes", "known_TS_genes",  "OMIM",
				"ClinGen_regions", "DGV", "DECIPHER", "gnomAD", "ClinVar_pathogenic", "ClinVar_benign", 
				"pathogenicity", "total_points"]
	
	for key in SUGGESTED_SCORES["DEL"].keys():
		key = key[0:2]
		if key not in header_list:
			header_list.append(key.split("-")[0])

	header_line = f"{separator}".join(raw_header_list + header_list)

	return header_line


def run_in_parallel(function, params_list, cores):
	"""Runs a function in parallel.
	Args:
		function: A function to run.
		params_list: A list of arguments for the function.
		cores: Number of threads, taken from the command line arguments.

	Returns:
		returncodes: A list of return codes.
	"""
	params_num = len(params_list)
	tasks = min([int(cores), params_num])
	results = [None] * params_num
	with Pool(tasks) as pool:
		# async_results = [pool.apply_async(function, params) for params in params_list]
		# for i in tqdm(range(params_num), desc='Processed CNVs'):
		# 	results[i] = async_results[i].get()
		for result in list(tqdm(pool.imap(function, params_list), total=len(params_list), desc='Processed CNVs')):
			results.append(result)
	
	return results
		

def set_relation(dict, lst):

	keys = list(dict.keys())
	
	if not keys:
		return False
	
	for key in keys:
		if set(lst).issubset(set(eval(key))):
			return [key, "sub"]
		elif set(lst).issuperset(set(eval(key))):
			return [key, "super"]
		else:
			return False


def count_procoding_genes(gene_list):
	count = 0
	gene_family = {}
	for gene_info in gene_list:
		gene_name, gene_group_id, phenotype = gene_info
		if not gene_group_id or phenotype:
			count += 1
			continue
		
		result = set_relation(gene_family, gene_group_id)
		if result:
			key, relation = result
			if relation == "sub":
				gene_family[key].append(gene_name)
			else:
				lst = gene_family.pop(key)
				lst.append(gene_name)
				gene_family[str(gene_group_id)] = lst
				
		else:
			gene_family[str(gene_group_id)] = list()
			gene_family[str(gene_group_id)].append(gene_name)	
	
	return len(gene_family) + count


def PVS1_Strength(cnv, tx_refseq_name):
	
	pvs1_del = {
		"VeryStrong":1, "Strong": 2, "Moderate": 3, "Supporting": 4, "Unmet": 5, "Unset": 5
	}
	pvs1_dup = {
		"VeryStrong":1, "Strong": 2, "Unmet": 3, "Unset": 3
	}

	cnv_info = CNVRecord(cnv.chr, cnv.start, cnv.end, cnv.type)
	if cnv.build == "hg38":
		transcripts = transcripts_hg38
	else:
		transcripts = transcripts_hg19
		
	tx_refseq_name_lst = tx_refseq_name.split("|")
	for tx_refseq in tx_refseq_name_lst:
		tx = get_transcript(tx_refseq, transcripts)
		if tx:
			break
	if not tx:
		print(tx_refseq_name)
		if cnv.type == "DEL":
			return "Unset",  5
		else:
			return "Unset",  3
		
	pvs1 = PVS1CNV(cnv_info, None, tx, cnv.build)
	if cnv.type == "DEL":
		strength = str(pvs1.verify_DEL().name)
		return strength, pvs1_del[strength]
	else:
		strength = str(pvs1.verify_DUP().name)
		return strength, pvs1_dup[strength]

def vcf2bed(vcf, bed):

	out = open(bed, "w")
	print("\t".join(["#chr", "start", "end", "type", "info"]), file = out)
	for line in open(vcf, "r"):
		if line.startswith("#"):
			continue
		line = line.rstrip().split("\t")

		if ";SVTYPE=" in line[7]:
			sv_type = line[7].split(";SVTYPE=")[1].split(";")[0]
		else:
			sv_type = line[4][1:-1]
			
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


def CNV_record_normalize(record):
	chr = record.split(":")[0]
	st, fl, type = record.split(":")[-1].split("-")
	type = type.upper()
	if not chr.startswith("chr"):
		chr = "chr" + chr
	cnv_list = [chr, st ,fl, type]
	return cnv_list


def CNV_record_build_convert(CNV, tmp_dir, chain_file):
	cnv_record = "-".join(CNV)
	src_bed = os.path.join(tmp_dir, f"{cnv_record}.bed")
	print("\t".join(CNV), file = open(src_bed, "w"))
	converted_bed = os.path.join(tmp_dir, os.path.basename(src_bed).replace(".bed", ".hg38.converted.bed"))
	unmap_bed = converted_bed + ".unmap"
	map_ratio = 0.95
	os.system(f"source activate interpret_env; CrossMap region {chain_file} {src_bed} {converted_bed} -r {map_ratio}")
	if not os.path.getsize(unmap_bed):
		print(f"Conversion success, map_ratio threshold is {map_ratio}")
		return converted_bed
	else:
		raise Exception("Conversion failed due to unmet map_ratio or unsuccessful mapping, please check the input CNV record")


def file_build_convert(file, tmp_dir, chain_file):
	converted_bed = os.path.join(tmp_dir, os.path.basename(file).replace(".bed", ".hg38.converted.bed"))
	unmap_bed = converted_bed + ".unmap"
	map_ratio = 0.95
	os.system(f"source activate interpret_env; CrossMap region {chain_file} {file} {converted_bed} -r {map_ratio}")
	if os.path.getsize(converted_bed):
		print(f"Conversion success, map_ratio threshold is {map_ratio}\n")
		print(f"The unmap CNV records is saved in: {unmap_bed}\n")
		return converted_bed
	else:
		raise Exception("Conversion failed due to unmet map_ratio or unsuccessful mapping, please check the input CNV record")
