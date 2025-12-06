import sys


infile, sample_name, outfile = sys.argv[1:4]
out = open(outfile, "w")


out.write('##fileformat=VCFv4.1\n')
out.write('##source=CNVpytor\n')
out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
out.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
out.write('##INFO=<ID=pytorRD,Number=1,Type=Float,Description="Normalized RD">\n')
out.write('##INFO=<ID=pytorP1,Number=1,Type=Float,Description="e-val by t-test">\n')
out.write('##INFO=<ID=pytorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">\n')
out.write('##INFO=<ID=pytorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">\n')
out.write('##INFO=<ID=pytorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">\n')
out.write('##INFO=<ID=pytorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">\n')
out.write('##INFO=<ID=pytorPE,Number=1,Type=Integer,Description="Number of paired-ends support the event">\n')
out.write('##INFO=<ID=pN,Number=1,Type=Float,Description="Fraction of reference genome gaps (Ns) in call region">\n')
out.write('##INFO=<ID=dG,Number=1,Type=Float,Description="distance from closest large (>100bp) gap in reference genome">\n')
out.write('##ALT=<ID=DEL,Description="Deletion">\n')
out.write('##ALT=<ID=DUP,Description="Duplication">\n')
out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
out.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
out.write('##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-ends that support the event">\n')
out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")


for line in open(infile, "r"):
    lst = line.rstrip().split("\t")
    type, id, length, rd, p1, p2, p3, p4, q0, pN, dG = lst
    chr = id.split(":")[0]
    st = id.split(":")[-1].split("-")[0]
    fl = id.split(":")[-1].split("-")[-1]
    type = type[0:3].upper()
    GT = ""
    if type == "DEL" and float(rd) < 0.25:
        GT = "1/1:0"
    elif type == "DEL" and float(rd) >= 0.25:
        GT = "0/1:1"
    elif type == "DUP" and float(rd) <= 1.75:
        GT = "0/1:3"
    elif type == "DUP" and float(rd) > 1.75:
        cn = round(2 * float(rd)) 
        GT = f"./1:{cn}"
    else:
        GT = "./.:."

    info = f"END={fl};SVTYPE={type};SVLEN={length};IMPRECISE;pytorRD={rd};pytorP1={p1};pytorP2={p2};pytorP3={p3};pytorP4={p4};pytorQ0={q0};pN={pN};dG={dG}"
    type = f"<{type}>"
    print("\t".join([chr, st, ".", "N", type, ".", "PASS", info, "GT:CN", GT]), file = out)

out.close()
