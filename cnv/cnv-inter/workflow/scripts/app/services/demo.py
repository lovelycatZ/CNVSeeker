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

print(retrieval_records_1(["7p15.2", "7p15.1", "7p14.3", "7p14.2", "7p14.1"], "DUP"))
print(retrieval_records_2(["7p15.2", "7p15.1", "7p14.3", "7p14.2", "7p14.1"], "chr7", "DUP"))
