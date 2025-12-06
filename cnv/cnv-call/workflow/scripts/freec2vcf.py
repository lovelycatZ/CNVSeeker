import sys


infile, sample_name, chrom_style, outfile = sys.argv[1:5]

with open(infile, "r") as f, open(outfile, "w") as out:
    next(f)
    out.write('##fileformat=VCFv4.1\n')
    out.write('##source=ControlFREEC\n')
    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    out.write('##INFO=<ID=W,Number=1,Type=Float,Description="Wilcoxon Rank Sum Test Pvalue">\n')
    out.write('##INFO=<ID=KS,Number=1,Type=Float,Description="Kolmogorov Smirnov Pvalue">\n')
    out.write('##ALT=<ID=DEL,Description="Deletion">\n')
    out.write('##ALT=<ID=DUP,Description="Duplication">\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Integer copy number">\n')
    out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
    for line in f:
        lst = line.rstrip().split("\t")
        chr, st, fl, cn, cnv_type, W, KS = lst
        st = int(st) + 1
        cn = int(cn)
        if chrom_style == "s" and chr.startswith("chr"):
            chr = chr[3:]
        elif chrom_style == "l" and not chr.startswith("chr"):
            chr = "chr" + chr
        else:
            pass

        if cn == 0:
            gt = "1/1"
        elif cn in [1, 3]:
            gt = "0/1"
        elif cn > 3:
            gt = "./."
        else:
            continue
        length = int(fl) - int(st) + 1
        cnv_type = "DEL" if cnv_type == "loss" else "DUP"
        info = f"SVTYPE={cnv_type};SVLEN={length};END={fl};W={W};KS={KS}"
        print("\t".join([chr, str(st), ".", "N", f"<{cnv_type}>", ".", "PASS", info, "GT:CN", f"{gt}:{cn}"]), file = out)
