import sys


bed, sample_name, vcf = sys.argv[1:4]
with open(bed, "r") as f, open(vcf, "w") as out:
    out.write('##fileformat=VCFv4.1\n')
    out.write('##source=ECOLE\n')
    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    out.write('##ALT=<ID=DEL,Description="Deletion">\n')
    out.write('##ALT=<ID=DUP,Description="Duplication">\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
    for line in f:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        _, chr, st, fl, cnv_type = line
        if cnv_type not in ["DUP", "DEL"]:
            continue
        sv_len = int(fl) - int(st)
        sv_info = f"END={fl};SVTYPE={cnv_type};SVLEN={sv_len}"
        print("\t".join([chr, str(int(st) + 1), ".", "N", f"<{cnv_type}>", ".", "PASS", sv_info, "GT", "./."]), file = out)
        