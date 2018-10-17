#!/usr/bin/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

import os
import sys
from docopt import docopt
sdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, sdir)
from produce_fa import produce_fasta
from art_produce_fq import produce_fq
from replacement_bp import replacement_bp
from merge_file import merge_file
from mapping import mapping
from rep_freebayes import freebayes_testing
from rep_varscan2 import varscan_testing
from rep_GATK4 import GATK_testing
from rep_lumpy import lumpy_testing
from rep_delly import delly_testing
from rep_gridss import gridss_testing
from rep_pindel import pindel_testing

def main_replacement(reads_length, multiple_count, single_count, bp, ti, multi, freq, tools, result_info):

    raw_fasta_prefix = "raw_replacement_" + str(ti + 1) + "_" + str(bp) + "bp"
    raw_fasta = produce_fasta(raw_fasta_prefix)
    multiple_output = "/".join(raw_fasta.split("/")[:2]) + "/multiple" + '_' + str(ti+1) + '_rep_' + str(bp) + "bp" + "_"
    single_output = "/".join(raw_fasta.split("/")[:2]) + "/single" + '_' + str(ti+1) + '_rep_' + str(bp) + "bp" + "_"
    produce_fq(raw_fasta, reads_length, multiple_count, multiple_output)
    info_record = replacement_bp(bp, raw_fasta)
    new_fasta = info_record['fasta']
    info_record['raw_fasta'] = raw_fasta
    info_record['multiple'] = multi
    produce_fq(new_fasta, reads_length, single_count, single_output)
    files = [multiple_output, single_output]
    fq1, fq2 = merge_file(files)

    print(info_record)
    ### add tools to be tested here ###
    bam_out, flag_out, stats_out = mapping(raw_fasta, fq1, fq2)

    call_tools = {'gatk':GATK_testing, 'freebayes':freebayes_testing, 'varscan':varscan_testing, 'lumpy':lumpy_testing, 'delly':delly_testing, 'gridss':gridss_testing, 'pindel':pindel_testing, }

    #### if no tools input, execute all tools ####
    if not tools:
        sum_tools = list(call_tools.keys())
    else:
        sum_tools = tools.lower().split(",")

    for tool in sum_tools:
        result_info += call_tools[tool](raw_fasta, bam_out, info_record, freq)

    ### delete the redundant files ###
    bam_index = bam_out + ".bai"
    new_fa_json = ".".join(new_fasta.split(".")[:-1]) + ".json"
    raw_fa_bwt = raw_fasta + ".bwt"
    raw_fa_pac = raw_fasta + ".pac"
    raw_fa_ann = raw_fasta + ".ann"
    raw_fa_amb = raw_fasta + ".amb"
    raw_fa_sa = raw_fasta + ".sa"
    raw_fa_dict = ".".join(raw_fasta.split(".")[:-1]) + ".dict"
    raw_fa_fai = raw_fasta + ".fai"
    del_files = [fq1,fq2,bam_out,flag_out,stats_out,bam_index,new_fa_json,raw_fa_bwt,raw_fa_pac,raw_fa_ann,raw_fa_amb,raw_fa_sa,raw_fa_dict,raw_fa_fai]
    for de in del_files:
        if os.path.exists(de):
            os.remove(de)

    return result_info

def parse_parameters(arguments):

    reads_length = arguments['--length']
    multiple_count = int(arguments['--count']) * (1 - float(arguments['MULTIPLE']))
    single_count = int(arguments['--count']) * float(arguments['MULTIPLE'])
    bp = arguments['--basepair']
    repeat_time = arguments['--times']
    multi = arguments['MULTIPLE']
    freq = arguments['--min_freq']
    tools = arguments['--tools']

    ###### use the function to test tools ######
    result_info = ""
    title = "Type" + "\t" + "Var_Type" + "\t" + "Interval" + "\t" + "Raw" + "\t" + "New" + "\t" + "Fasta" + "\t" + "Multiple" + "\t" + "Tool" + "\t" + "Tool_POS" + "\t" + "Tool_REF" + "\t" + "Tool_ALT" + "\t" + "Tool_Type" + "\t" + "Tool_ratio" + "\n"

    for ti in range(int(repeat_time)):
        result_info += main_replacement(reads_length, multiple_count, single_count, bp, ti, multi, freq, tools, result_info)

    rep_res = "./data/simulate" + "_rep_" + bp + "_tools.txt"
    with open(rep_res, 'w') as fw:
        fw.write(title)
        fw.write(result_info)

if __name__ == "__main__":
    usage = """
    Usage:
        simulate_replacement.py [-l=150] [-c=100000] [-b=1] [-t=1] [-v=0.2] [--tools=<arg>] MULTIPLE

    Testing different tools on different raw-fasta-based replacement variation

    Arguments:
        MULTIPLE        the proportion of replaced-variation fastas in all fastas

    Options:
        -h --help
        -l,--length=150         reads length of simulated fasta [default: 150]
        -c,--count=100000       number of reads/read pairs [default: 100000]
        -b,--basepair=1         the length of replaced bases [default: 1]
        -t,--times=1            the repeat time [default: 1]
        -v,--min_freq=0.2       Minimum variant allele frequency threshold [default: 0.2]
        --tools=<arg>           identify the tested tool, available values are: Freebayes,Varscan, Delly, Gridss, Pindel, GATK and Lumpy
    """

    arguments = docopt(usage)
    parse_parameters(arguments)
