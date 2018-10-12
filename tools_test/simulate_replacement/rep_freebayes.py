#!/bin/usr/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

try:
    from config import freebayes
except:
    freebayes = "/Users/xu/anaconda3/bin/freebayes"

import os
import sys
import subprocess

def merge_files(info_record, new):

    with open(new, 'r') as fn:

        ori_loc = info_record['info'][1].split('-')[0]
        Type = "Replacement"
        Var_Type = info_record['info'][0]
        Interval = info_record['info'][1]
        Raw = info_record['info'][2]
        New = info_record['info'][3]
        fasta = info_record['raw_fasta']
        multiple = info_record['multiple']
        ori_info = Type + "\t" + Var_Type + "\t" + Interval + "\t" + Raw + "\t" + New + "\t" + fasta + "\t" + str(multiple) + "\t" + "Freebayes"

        mes = ""
        for fn_line in fn:
            if "#" not in fn_line:
                pos = fn_line.strip().split()[1]
                ref = fn_line.strip().split()[3]
                alt = fn_line.strip().split()[4]
                var_type = fn_line.strip().split()[7].split(";")[40]
                ref_reads = fn_line.strip().split()[9].split(":")[3]
                alt_reads = fn_line.strip().split()[9].split(":")[5]
                ratio = int(alt_reads)/(int(ref_reads) + int(alt_reads))
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"
        
        if not mes:
            mes = ori_info + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\n"

    return mes

def freebayes_testing(fa, bam, info_record, *args):

    vcf = ".".join(bam.split(".")[:-1]) + "_freebayes.vcf"
    cmd = "%s --fasta-reference %s %s > %s 2>/dev/null" % (freebayes, fa, bam, vcf)
    subprocess.call(cmd, shell=True)
    mes = merge_files(info_record, vcf)

    return mes
