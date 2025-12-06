import os
import sys
import subprocess



def resolve_vcf_list(vcf_list):
    vcf_dict = {}
    with open(vcf_list, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "cn.mops" in line:
                tool = "cn.mops"
            else:
                tool = os.path.basename(line).split(".")[1]
            vcf_dict[tool] = line
    return vcf_dict


# "(CNVpytor,Manta,Wham)(2) + (lumpy)(1)"

def resolve_combo(combo):
    combo_dict = {}
    for part in combo.split("+"):
        part = part.strip()
        if not part:
            continue
        tools_section, support_section = part.split(")(")
        tools_section = tools_section.replace("(", "").replace(")", "")
        support_section = support_section.replace("(", "").replace(")", "")
        tools = tools_section.split(",")
        support = int(support_section)
        combo_dict[part] = {"tools": tools, "support": support}

    return combo_dict


def make_vcf_list_file(vcf_list, vcf_list_file):
    with open(vcf_list_file, "w") as f:
        for vcf_path in vcf_list:
            f.write(vcf_path + "\n")


def run_survivor_merge(vcf_list, output_vcf, min_support, max_dist_ratio = 0.5):
    """调用 SURVIVOR 进行合并"""
    cmd = [
        "SURVIVOR", "merge", vcf_list,
        str(max_dist_ratio), str(min_support), "1", "1", "0", "50", output_vcf
    ]
    subprocess.run(cmd, check=True)


def main():
    vcf_list = sys.argv[1]
    output_vcf = sys.argv[2]
    tmp_dir = sys.argv[3]
    combo = sys.argv[4]

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    vcf_dict = resolve_vcf_list(vcf_list)
    combo_dict = resolve_combo(combo)

    # 对于每套组合，生成中间结果
    intermediate_vcfs = []
    for combo_name, info in combo_dict.items():
        tools = info["tools"]
        min_support = info["support"]
        vcf_list = [vcf_dict[tool] for tool in tools]
        vcf_list_file = os.path.join(tmp_dir, f"vcf_list_{combo_name}.txt")
        out_file = os.path.join(tmp_dir, f"merged_{combo_name}.vcf")
        make_vcf_list_file(vcf_list, vcf_list_file)
        run_survivor_merge(vcf_list_file, out_file, min_support)
        intermediate_vcfs.append(out_file)

    # 最后将所有中间结果合并
    final_vcf_list = os.path.join(tmp_dir, "final_vcf_list.txt")
    with open(final_vcf_list, "w") as f:
        for vcf in intermediate_vcfs:
            f.write(vcf + "\n")
    run_survivor_merge(final_vcf_list, output_vcf, 1)


if __name__ == "__main__":
    main()
