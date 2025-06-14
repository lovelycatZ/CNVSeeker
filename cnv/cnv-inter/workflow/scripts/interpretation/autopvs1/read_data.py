#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27
# modified: 2021/6/29

import os
import configparser
from pyfaidx import Fasta
from .pyhgvs.utils import read_transcripts
from .utils import read_pvs1_levels, create_bed_dict

BinPath = os.path.split(os.path.realpath(__file__))[0]

config = configparser.ConfigParser()
config.read(BinPath+'/config.ini')


for top in config:
    for key in config[top]:
        if not config[top][key].startswith('/'):
            config[top][key] = os.path.join(BinPath, config[top][key])
    

pvs1_levels = read_pvs1_levels(config['DEFAULT']['pvs1levels'])


# genome_hg19 = Fasta(config['HG19']['genome'])
# genome_hg38 = Fasta(config['HG38']['genome'])

transcripts_hg19 = read_transcripts(open(config['HG19']['transcript']))
transcripts_hg38 = read_transcripts(open(config['HG38']['transcript']))

domain_hg19 = create_bed_dict(config['HG19']['domain'])
domain_hg38 = create_bed_dict(config['HG38']['domain'])

hotspot_hg19 = create_bed_dict(config['HG19']['hotspot'])
hotspot_hg38 = create_bed_dict(config['HG38']['hotspot'])

curated_region_hg19 = create_bed_dict(config['HG19']['curated_region'])
curated_region_hg38 = create_bed_dict(config['HG38']['curated_region'])

exon_lof_popmax_hg19 = create_bed_dict(config['HG19']['exon_lof_popmax'])
exon_lof_popmax_hg38 = create_bed_dict(config['HG38']['exon_lof_popmax'])

# pathogenic_hg19 = read_pathogenic_site(config['HG19']['pathogenic_site'])
# pathogenic_hg38 = read_pathogenic_site(config['HG38']['pathogenic_site'])
