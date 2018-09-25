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

def main_replacement(reads_length, multiple_count, single_count, bp, repeat_time, multi, freq):

    result_info = ""
    title = "Type" + "\t" + "Var_Type" + "\t" + "Interval" + "\t" + "Raw" + "\t" + "New" + "\t" + "Fasta" + "\t" + "Multiple" + "\t" + "Tool" + "\t" + "Tool_POS" + "\t" + "Tool_REF" + "\t" + "Tool_ALT" + "\t" + "Tool_Type" + "\t" + "Tool_ratio" + "\n"

    result_info += title

    for ti in range(int(repeat_time)):
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

        ### add tools to be tested here ###
        bam_out, flag_out, stats_out = mapping(raw_fasta, fq1, fq2)

        ##### test tools of FreeBayes #####
        freebayes_out = freebayes_testing(raw_fasta, bam_out, info_record)
        result_info += freebayes_out

        ##### test tools of VarScan2 ######
        varscan_out = varscan_testing(raw_fasta, bam_out, info_record, freq)
        result_info += varscan_out

        ##### test tools of GATK4 ######
        gatk_out = GATK_testing(raw_fasta, bam_out, info_record)
        result_info += gatk_out

        ##### test tools of lumpy ######
        lumpy_out = lumpy_testing(bam_out, info_record)
        result_info += lumpy_out
        ###################################

    rep_res = bam_out.split("_")[0] + "_rep_" + bp + "_tools.txt"
    with open(rep_res, 'w') as fw:
        fw.write(result_info)

    return rep_res

def parse_parameters(arguments):

    reads_length = arguments['--length']
    multiple_count = int(arguments['--count']) * (1 - float(arguments['MULTIPLE']))
    single_count = int(arguments['--count']) * float(arguments['MULTIPLE'])
    bp = arguments['--basepair']
    repeat_time = arguments['--times']
    multi = arguments['MULTIPLE']
    freq = arguments['--min_freq']

    ###### use the function to test tools ######
    main_replacement(reads_length, multiple_count, single_count, bp, repeat_time, multi, freq)

if __name__ == "__main__":
    usage = """
    Usage:
        simulate_replacement.py [-l=150] [-c=100000] [-b=1] [-t=1] [-v=0.2] MULTIPLE

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
    """

    arguments = docopt(usage)
    parse_parameters(arguments)
