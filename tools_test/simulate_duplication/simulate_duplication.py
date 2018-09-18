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
from duplication_bp import duplication_bp
from merge_file import merge_file
from mapping import mapping
from dup_freebayes import freebayes_testing

def main_duplication(reads_length, multiple_count, single_count, bp, repeat_time, copy_number, multi):

    result_info = ""
    title = "Type" + "\t" + "Var_Type" + "\t" + "Interval" + "\t" + "Raw" + "\t" + "New" + "\t" + "Fasta" + "\t" + "Multiple" + "\t" + "Tool" + "\t" + "Tool_POS" + "\t" + "Tool_REF" + "\t" + "Tool_ALT" + "\t" + "Tool_Type" + "\t" + "Tool_ratio" + "\n"

    result_info += title

    for ti in range(int(repeat_time)):
        raw_fasta_prefix = "raw_duplication_" + str(ti + 1) + "_"+ str(bp) + "bp"
        raw_fasta = produce_fasta(raw_fasta_prefix)
        multiple_output = "/".join(raw_fasta.split("/")[:2]) + "/multiple" + '_' + str(ti+1) + '_dup_' + str(bp) + "bp" + "_"
        single_output = "/".join(raw_fasta.split("/")[:2]) + "/single" + '_' + str(ti+1) + '_dup_' + str(bp) + "bp" + "_"
        produce_fq(raw_fasta, reads_length, multiple_count, multiple_output)
        info_record = duplication_bp(bp, raw_fasta, copy_number)
        new_fasta = info_record['fasta']
        info_record['raw_fasta'] = raw_fasta
        info_record['multiple'] = multi
        produce_fq(new_fasta, reads_length, single_count, single_output)
        files = [multiple_output, single_output]
        fq1, fq2 = merge_file(files)
        ### add tools to be tested here ###

        bam_out, flag_out, stats_out = mapping(raw_fasta, fq1, fq2)
        freebayes_out = freebayes_testing(raw_fasta, bam_out, info_record)
        result_info += freebayes_out

        ###################################

    dup_res = bam_out.split("_")[0] + "_dup_" + bp + "_tools.txt"
    with open(dup_res, 'w') as fw:
        fw.write(result_info)

    return dup_res

def parse_parameters(arguments):

    reads_length = arguments['--length']
    multiple_count = int(arguments['--count']) * (1 - float(arguments['MULTIPLE']))
    single_count = int(arguments['--count']) * float(arguments['MULTIPLE'])
    bp = arguments['--basepair']
    repeat_time = arguments['--times']
    copy_number = arguments['--copys']
    multi = arguments['MULTIPLE']

    ###### use the function to test tools ######
    main_duplication(reads_length, multiple_count, single_count, bp, repeat_time, copy_number, multi)

if __name__ == "__main__":
    usage = """
    Usage:
        simulate_duplication.py [-l=150] [-c=100000] [-b=1] [-p=1] [-t=1] MULTIPLE

    Testing different tool on different raw-fasta-based duplication variation

    Arguments:
        MULTIPLE        the proportion of duplicated-variation fastas in all fastas

    Options:
        -h --help
        -l,--length=150         reads length of simulated fasta [default: 150]
        -c,--count=100000       number of reads/read pairs [default: 100000]
        -b,--basepair=1         the length of duplicated bases [default: 1]
        -p,--copys=1            the copy number [default: 1] 
        -t,--times=1            the repeat time [default: 1]
    """

    arguments = docopt(usage)
    parse_parameters(arguments)
