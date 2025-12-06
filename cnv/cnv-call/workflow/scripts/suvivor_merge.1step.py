import os
import sys
import subprocess


COMBO_DICT = {
	# "HD_WGS":{
	# 	"DEL":{
	# 		"SS": "(Wham,delly)(2)+(lumpy)(1)", # yes
	# 		"S": "(delly,lumpy)(2)", # yes
	# 		"M": "(CNVpytor,cn.mops,delly,lumpy)(2)", # yes
	# 		# "L": "(CNVpytor,ControlFREEC,cn.mops)(2)" # yes
	# 		"L": "(CNVpytor,delly,lumpy)(3)+(ControlFREEC,cn.mops)(2)"
	# 	},
	# 	"DUP":{
	# 		"SS": "(Wham)(1)+(delly,lumpy)(2)", # yes
	# 		"S": "(Wham,lumpy)(2)+(delly)(1)", # yes
	# 		"M": "(CNVpytor,ControlFREEC)(1)+(Wham,cn.mops,delly,lumpy)(2)", # yes
	# 		"L": "(CNVpytor)(1)+(ControlFREEC,cn.mops)(2)", # yes
	# 	}
	# },
	"HD_WGS":{
		"DEL":{
			"SS": "(Wham,delly)(2)+(lumpy)(1)", # yes
			"S": "(delly,lumpy)(2)", # yes
			"M": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)", # yes
            "L": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)"
		},
		"DUP":{
			"SS": "(Wham)(1)+(delly,lumpy)(2)", # yes
			"S": "(Wham,lumpy)(2)+(delly)(1)", # yes
			"M": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)", # yes
            "L": "(CNVpytor,ControlFREEC,Wham,cn.mops,delly,lumpy)(3)"
		}
	},
	"LD_WGS":{
		"DEL":{
			"SS": "(delly,gridss,lumpy)(1)", # yes
			"S": "(delly,gridss,lumpy)(1)", # yes
			"M": "(CNVpytor,cn.mops,delly,lumpy)(2)+(gridss)(1)", # yes
			"L": "(CNVpytor,delly)(2)+(ControlFREEC,cn.mops)(2)" # yes
		},
		"DUP":{
			"SS": "(Manta,Wham,gridss)(1)", # yes
			"S": "(Manta,gridss)(1)", # yes
			"M": "(CNVpytor,Manta,Wham)(1)+(ControlFREEC,gridss)(1)", # yes
			"L": "(CNVpytor)(1)+(ControlFREEC,Manta,Wham,cn.mops)(3)", # yes
		}
	},
	"WES":{
		"DEL":{
			"S": "(ECOLE,ExomeDepth,gatk4,xhmm)(2)+(cn.mops)(1)", # yes
			"M": "(ECOLE,ExomeDepth,cn.mops,gatk4,xhmm)(1)", # yes
			"L": "(ECOLE,gatk4)(1)+(ExomeDepth,xhmm)(1)" # yes
		},
		"DUP":{
			"S": "(ECOLE)(1)+(ExomeDepth)(1)", # yes
			"M": "(ECOLE)(1)+(ExomeDepth,gatk4,xhmm)(2)", # yes
			"L": "(ECOLE,ExomeDepth,cn.mops,xhmm)(1)", # yes
		}
	}
}


SIZE_RANGE_DICT = {
    "1": {
        "S": (50, 10000),
        "M": (10000, 100000),
        "L": (100000, 10000000)
    },
    "2": {
        "SS": (50, 500),
        "S": (500, 5000),
        "M": (5000, 100000),
        "L": (100000, 10000000)
    }
}


# def resolve_vcf_list(vcf_list):
#     vcf_dict = {}
#     with open(vcf_list, "r") as f:
#         for line in f:
#             line = line.strip()
#             if not line:
#                 continue
#             if "cn.mops" in line:
#                 tool = "cn.mops"
#             else:
#                 tool = os.path.basename(line).split(".")[1]
#             vcf_dict[tool] = line
#     return vcf_dict


def resolve_vcf_list(vcf_list):
    vcf_dict = {}
    with open(vcf_list, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            tool = os.path.basename(line).split(".vcf")[0].split("_")[-1]
            vcf_dict[tool] = line
    return vcf_dict


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


def get_standard_vcf(input_vcf, output_vcf, tool, sample_id, cnv_type, min_len=50, max_len=10000000, apply_filter=True, sex='M'):
    script = "/gpfs/hpc/home/lijc/xiangxud/project/new_result/hd_wgs/cnvseeker/standardize_vcf.py"
    cmd = f"python {script} {input_vcf} {output_vcf} {tool} {sample_id} {cnv_type} {min_len} {max_len} {str(apply_filter)} {sex}"
    os.system(cmd)


def main():
    vcf_list = sys.argv[1]
    final_merged_vcf = sys.argv[2]
    result_dir = sys.argv[3]
    data_type = sys.argv[4]
    sample_id = sys.argv[5]

    standard_vcf_dir = os.path.join(result_dir, "standard_vcfs")
    if not os.path.exists(standard_vcf_dir):
        os.makedirs(standard_vcf_dir, exist_ok=True)

    if data_type == 'WES':
        SIZE_RANGE = SIZE_RANGE_DICT['1']
    else:
        SIZE_RANGE = SIZE_RANGE_DICT['2']

    # 解析输入的vcf列表
    tool_vcf_dict = resolve_vcf_list(vcf_list)
    print(tool_vcf_dict)
    final_vcf_list = []
    for cnv_type in COMBO_DICT[data_type]:
        for bin_tag, combo_name in COMBO_DICT[data_type][cnv_type].items():
            # 解析组合方案
            combo_dict = resolve_combo(combo_name)

            # 创建临时目录
            merge_tmp_dir = os.path.join(result_dir, f"merge_tmp.{cnv_type}.{bin_tag}")
            if not os.path.exists(merge_tmp_dir):
                os.makedirs(merge_tmp_dir)

            # 对于每套组合，生成中间结果 
            intermediate_vcfs = []
            for combo_name, info in combo_dict.items():
                tools = info["tools"]
                min_support = info["support"]
                vcf_list = []
                for tool in tools:
                    tool_input_vcf = tool_vcf_dict.get(tool, "")
                    if not tool_input_vcf:
                        exit(f"Error: VCF for tool {tool} not found in the provided VCF list.")

                    tool_standard_vcf = os.path.join(standard_vcf_dir, f"{sample_id}.{tool}.{cnv_type}.{bin_tag}.vcf")
                    size_range = SIZE_RANGE[bin_tag]
                    min_len, max_len = size_range              
                    get_standard_vcf(tool_input_vcf, tool_standard_vcf, tool, sample_id, cnv_type, min_len, max_len, True, 'M')
                    vcf_list.append(tool_standard_vcf)

                tmp_vcf_list = os.path.join(merge_tmp_dir, f"vcf_list_{combo_name}.txt")
                tmp_vcf = os.path.join(merge_tmp_dir, f"merged_{combo_name}.vcf")
                make_vcf_list_file(vcf_list, tmp_vcf_list)
                run_survivor_merge(tmp_vcf_list, tmp_vcf, min_support)
                intermediate_vcfs.append(tmp_vcf)

            # 最后将所有中间结果合并
            tmp_final_vcf_list = os.path.join(merge_tmp_dir, "final.vcf_list.txt")
            tmp_final_vcf = os.path.join(merge_tmp_dir, f"{sample_id}.{cnv_type}.{bin_tag}.merged.vcf")
            with open(tmp_final_vcf_list, "w") as f:
                for vcf in intermediate_vcfs:
                    f.write(vcf + "\n")
            run_survivor_merge(tmp_final_vcf_list, tmp_final_vcf, 1)
            final_vcf_list.append(tmp_final_vcf)

    # 最终将所有类型和大小段的结果合并
    final_vcf_list_file = os.path.join(result_dir, f"{sample_id}.vcf_list.txt")
    make_vcf_list_file(final_vcf_list, final_vcf_list_file)
    run_survivor_merge(final_vcf_list_file, final_merged_vcf, 1)



if __name__ == "__main__":
    main()
