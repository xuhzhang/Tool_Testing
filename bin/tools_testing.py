#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from docopt import docopt
from concurrent.futures import ProcessPoolExecutor
sdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../tools_test/")
sys.path.insert(0, sdir)
from simulate_deletion.simulate_deletion import main_deletion
from simulate_duplication.simulate_duplication import main_duplication
from simulate_insertion.simulate_insertion import main_insertion
from simulate_inversion.simulate_inversion import main_inversion
from simulate_replacement.simulate_replacement import main_replacement

def call_functions(arguments):

    reads_length = arguments['--length']
    reads_count = arguments['--count']
    bp = arguments['--basepair']
    repeat_time = arguments['--times']
    copy_number = arguments['--copys']
    freq = arguments['--min_freq']
    tools = arguments['--tools']
    worker_num = int(arguments['--thread'])

    sum_info = []
    work_pool = ProcessPoolExecutor(int(worker_num))

    if arguments['--dup']:
        dup_ratio = arguments['--dup']
        dup_multiple_count = int(reads_count) * (1 - float(dup_ratio))
        dup_single_count = int(reads_count) * float(dup_ratio)

        dup_in_info = ""
        for ti in range(int(repeat_time)):
            dup_info = work_pool.submit(main_duplication, reads_length, dup_multiple_count, dup_single_count, bp, ti, copy_number, dup_ratio, freq, tools, dup_in_info)
            sum_info.append(dup_info)

    if arguments['--del']:
        del_ratio = arguments['--del']
        del_multiple_count = int(reads_count) * (1 - float(del_ratio))
        del_single_count = int(reads_count) * float(del_ratio)

        del_in_info = ""
        for ti in range(int(repeat_time)):
            del_info = work_pool.submit(main_deletion, reads_length, del_multiple_count, del_single_count, bp, ti, del_ratio, freq, tools, del_in_info)
            sum_info.append(del_info)

    if arguments['--ins']:
        ins_ratio = arguments['--ins']
        ins_multiple_count = int(reads_count) * (1 - float(ins_ratio))
        ins_single_count = int(reads_count) * float(ins_ratio)

        ins_in_info = ""
        for ti in range(int(repeat_time)):
            ins_info = work_pool.submit(main_insertion, reads_length, ins_multiple_count, ins_single_count, bp, ti, ins_ratio, freq, tools, ins_in_info)
            sum_info.append(ins_info)

    if arguments['--inv']:
        inv_ratio = arguments['--inv']
        inv_multiple_count = int(reads_count) * (1 - float(inv_ratio))
        inv_single_count = int(reads_count) * float(inv_ratio)

        inv_in_info = ""
        for ti in range(int(repeat_time)):
            inv_info = work_pool.submit(main_inversion, reads_length, inv_multiple_count, inv_single_count, bp, ti, inv_ratio, freq, tools, inv_in_info)
            sum_info.append(inv_info)

    if arguments['--rep']:
        rep_ratio = arguments['--rep']
        rep_multiple_count = int(reads_count) * (1 - float(rep_ratio))
        rep_single_count = int(reads_count) * float(rep_ratio)

        rep_in_info = ""
        for ti in range(int(repeat_time)):
            rep_info = work_pool.submit(main_replacement, reads_length, rep_multiple_count, rep_single_count, bp, ti, rep_ratio, freq, tools, rep_in_info)
            sum_info.append(rep_in_info)

    work_pool.shutdown(wait=True)
    return sum_info

def main_function(arguments):

    sum_info = call_functions(arguments)

    return sum_info


if __name__ == "__main__":
    usage = """
    Usage:
        tools_testing.py [-l=150] [-c=100000] [-b=1] [-t=1] [-p=1] [-d=1] [-v=0.2] [--tools=<arg>] [--dup <dup-ration>] [--del <del-ration>] [--ins <ins-ration>] [--inv <inv-ration>] [--rep <rep-ration>]

    Testing different tools on different raw-fasta-based variations

    Options:
        -h --help
        -l,--length=150             reads length of simulated fasta [default: 150]
        -c,--count=100000           number of reads/read pairs [default: 100000]
        -b,--basepair=1             the length of variated bases [default: 1]
        -p,--copys=1                the copy number variation of duplication [default: 1]
        -t,--times=1                the repeat time [default: 1]
        -d,--thread=1               the worker threads [default: 1]
        -v,--min_freq=0.2           Minimum variant allele frequency threshold of VarScan2 [default: 0.2]
        --tools=<arg>               identify the tested tool, available values are: Freebayes, Varscan, Delly, Gridss, GATK and Lumpy
        --dup <dup-ratio>           the proportion of duplicated-variation reads in all reads
        --del <del-ratio>           the propertion of deleted-variation reads in all reads
        --ins <ins-ratio>           the propertion of inserted-variation reads in all reads
        --inv <inv-ratio>           the propertion of inversed-variation reads in all reads
        --rep <rep-ratio>           the propertion of repeated-variation reads in all reads
    """

    arguments = docopt(usage)
    sum_info = main_function(arguments)

    total_file = "./data/total_info_aggregate.txt"
    title = "Type" + "\t" + "Var_Type" + "\t" + "Interval" + "\t" + "Raw" + "\t" + "New" + "\t" + "Fasta" + "\t" + "Multiple" + "\t" + "Tool" + "\t" + "Tool_POS" + "\t" + "Tool_REF" + "\t" + "Tool_ALT" + "\t" + "Tool_Type" + "\t" + "Tool_ratio" + "\n"
    with open(total_file, 'w') as fw:
        fw.write(title)
        for info in sum_info:
            fw.write(info.result())
