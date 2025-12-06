import pandas as pd
import pysam
from collections import namedtuple
import sys


# 定义结构体来存储gap区域信息
GapRegion = namedtuple('GapRegion', ['chrom', 'start', 'end'])
CHROMS = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]


def load_gap_regions(gap_file):
    """
    加载gap区域文件。
    """
    gap_regions = []
    with open(gap_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                chrom = fields[0]
                if chrom not in CHROMS:
                    continue
                start = int(fields[1])
                end = int(fields[2])
                gap_regions.append(GapRegion(chrom, start, end))
    return gap_regions


def overlap_length(region1, region2):
    """
    计算两个区域重叠的长度。
    """
    if region1.chrom != region2.chrom:
        return 0
    
    start = max(region1.start, region2.start)
    end = min(region1.end, region2.end)
    if start < end:
        return end - start
    
    return 0


def filter_vcf(input_vcf, gap_regions, output_vcf):
    """
    根据gap区域过滤VCF文件，生成新的VCF文件。
    """
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):  # VCF头部信息直接写入输出
                outfile.write(line)
                continue

            # 处理变异行
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]  # INFO字段，第8列
            end = int(info.split('END=')[1].split(';')[0])  # 获取END位置
            sv_len = end - pos + 1
            
            variant_start = pos
            variant_end = end  # DUP的结束位置

            # 计算DEL与gap区域的重叠长度
            total_overlap_length = 0
            for gap in gap_regions:
                    total_overlap_length += overlap_length(gap, GapRegion(chrom, variant_start, variant_end))

            # 如果DEL的重叠部分超过DEL长度的70%，则过滤掉该DEL
            if total_overlap_length > 0.7 * sv_len:
                continue

            # 如果不满足过滤条件，将该变异行写入输出文件
            outfile.write(line)


if __name__ == "__main__":
    # 输入文件路径和输出文件路径
    input_vcf, gap_file, output_vcf = sys.argv[1], sys.argv[2], sys.argv[3]

    # 加载gap区域
    gap_regions = load_gap_regions(gap_file)

    # 过滤VCF文件
    filter_vcf(input_vcf, gap_regions, output_vcf)
