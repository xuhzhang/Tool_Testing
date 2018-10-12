#!/bin/usr/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

try:
    from config import delly, bcftools
except:
    delly = "/home/syngen/Software/delly/delly_v0.7.8_linux_x86_64bit"
    bcftools = "/home/syngen/Software/bcftools-1.9/bcftools"

import os
import sys
import subprocess

def merge_files(info_record, new):

    with open(new, 'r') as fn:
        ori_loc = info_record['info'][1].split('-')[0]
        Type = "Duplication"
        Var_Type = info_record['info'][0]
        Interval = info_record['info'][1]
        Raw = info_record['info'][2]
        New = info_record['info'][3]
        fasta = info_record['raw_fasta']
        multiple = info_record['multiple']
        ori_info = Type + "\t" + Var_Type + "\t" + Interval + "\t" + Raw + "\t" + New + "\t" + fasta + "\t" + str(multiple) + "\t" + "Delly"

        mes = ""
        for fn_line in fn:
            if "#" not in fn_line:
                pos = fn_line.strip().split("\t")[1]
                ref = fn_line.strip().split("\t")[3]
                alt = "-"
                var_type = fn_line.strip().split("\t")[7].split(";")[1]
                ratio = "-"
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

        if not mes:
            mes = ori_info + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\n"

    return mes

def delly_testing(fa, bam, info_record, *args):

    bcf = ".".join(bam.split(".")[:-1]) + "_delly.bcf"
    cmd1 = "%s call -g %s -o %s %s > /dev/null 2>&1" % (delly, fa, bcf, bam)

    subprocess.call(cmd1, shell=True)

    vcf = ".".join(bam.split(".")[:-1]) + "_delly.vcf"
    cmd2 = "%s view %s > %s 2>/dev/null" % (bcftools, bcf, vcf)

    subprocess.call(cmd2, shell=True)

    mes = merge_files(info_record, vcf)

    bcf_index = bcf + ".csi"
    os.remove(bcf)
    os.remove(bcf_index)

    return mes
