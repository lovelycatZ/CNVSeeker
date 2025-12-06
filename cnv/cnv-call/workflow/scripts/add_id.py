import sys


def add_var_id(in_vcf, out_vcf):
	sv_type_count = {}
	out = open(out_vcf, "w")
	out.write('##fileformat=VCFv4.2\n')
	out.write('##source=CNVSeeker\n')
	out.write('##ALT=<ID=DEL,Description="Deletion">\n')
	out.write('##ALT=<ID=DUP,Description="Duplication">\n')
	out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
	out.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
	out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n')
	out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
	with open(in_vcf, "r") as f:
		for line in f:
			if line.startswith("#"):
				continue
			line = line.rstrip().split("\t")
			line = line[0:8]
			sv_type = line[7].split("SVTYPE=")[1].split(";")[0]
			if sv_type not in sv_type_count:
				sv_type_count[sv_type] = 0
			else:
				sv_type_count[sv_type] += 1
			line[2] = f"CNVSeeker.{sv_type}.{sv_type_count[sv_type]}"
			sv_len = int(line[7].split("SVLEN=")[1].split(";")[0])
			st = int(line[1])
			if sv_len < 0:
				sv_len = -sv_len
			fl = st + sv_len - 1
			line[7] = f"SVTYPE={sv_type};SVLEN={sv_len};END={fl}"
			line[3] = "N"
			line[5] = "."
			line[6] = "."
			print("\t".join(line), file = out)
	out.close()
	

if __name__ == "__main__":
	input_vcf = sys.argv[1]
	output_vcf = sys.argv[2]
	add_var_id(input_vcf, output_vcf)
