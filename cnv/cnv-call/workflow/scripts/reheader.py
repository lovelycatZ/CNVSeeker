import sys


org_header = sys.argv[1]
new_header = sys.argv[2]

with open(org_header, "r") as f, open(new_header, "w") as out:
    for line in f:
        if any(x in line for x in ["AH:", "PM:"]):
            line = line.strip().split("\t")
            new_line = [item for item in line if not item.startswith(("AH:", "PM:"))]
            print("\t".join(new_line), file = out)
        else:
            out.write(line)
