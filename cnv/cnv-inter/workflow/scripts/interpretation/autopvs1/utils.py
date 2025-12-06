#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/2/2 16:28

import sys
from collections import namedtuple


VCFRecord = namedtuple('VCFRecord', ("chrom", "pos", 'ref', 'alt'))


def get_transcript(trans_name, transcripts):
    """
    Get a transcript using its name or a gene name
    :param trans_name:
    :param transcripts:
    :return: transcript
    """
    transcript = transcripts.get(trans_name)
    if not transcript:
        transcript = transcripts.get(trans_name.split('.')[0])
    if not transcript:
        transcript = None
    return transcript


def create_two_dim_dict(thedict, key1, key2, val):
    """
    add two dimension dict
    :param thedict:
    :param key1:
    :param key2:
    :param val:
    :return: dict
    """
    if key1 in thedict:
        thedict[key1].update({key2: val})
    else:
        thedict.update({key1: {key2: val}})

        
def create_bed_dict(bed):
    """
    read bed3 / bed6 / bed12 format as a dict
    :param bed:
    :return: bed dict
    """
    bed_dict = dict()
    try:
        with open(bed)as bed:
            for line in bed:
                record = line.strip().split("\t")
                chrom = record[0]
                if len(record) >= 12:
                    block_count = record[9]
                    block_sizes = record[10].split(",")
                    block_starts = record[11].split(",")
                    for i in range(int(block_count)):
                        start = int(record[1]) + int(block_starts[i])
                        end = start + int(block_sizes[i])
                        key = record[3] + '|' + str(start) + '-' + str(end)
                        create_two_dim_dict(bed_dict, key, "chrom", chrom)
                        create_two_dim_dict(bed_dict, key, "start", int(start))
                        create_two_dim_dict(bed_dict, key, "end", int(end))
                else:
                    start = int(record[1])
                    end = int(record[2])
                    if len(record) > 3:
                        key = record[3]
                    else:
                        key = chrom + ":" + str(start) + "-" + str(end)
                    create_two_dim_dict(bed_dict, key, "chrom", chrom)
                    create_two_dim_dict(bed_dict, key, "start", int(start))
                    create_two_dim_dict(bed_dict, key, "end", int(end))
        return bed_dict
    except Exception as e:
        sys.stderr.write(str(e))


def contained_in_bed(bed_dict, chrom, start, end):
    """
    Determine whether a variant is contained in the bed region.
    :param bed_dict:
    :param chrom:
    :param start:
    :param end:
    :return: Boolean value
    """
    chrom = str(chrom)
    if "chr" not in chrom:
        chrom = "chr" + chrom
    # max(start1, start2) < min(end1, end2)
    for key in bed_dict:
        if bed_dict[key]["chrom"] == chrom and \
                max(bed_dict[key]["start"], start) < min(bed_dict[key]["end"], end):
            return True, key
    return False


def read_pvs1_levels(file):
    """
    :param file: PVS1 Levels
    :return: dict
    """
    _pvs1_levels = {}
    try:
        with open(file) as fh:
            for line in fh:
                record = line.strip().split("\t")
                gene = record[0]
                level = record[1]
                if gene not in _pvs1_levels:
                    _pvs1_levels[gene] = level
    except Exception as err:
        sys.stderr.write(err)

    return _pvs1_levels
