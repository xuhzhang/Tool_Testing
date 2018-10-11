#!/usr/bin/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

try:
    from config import pindel, pindel2vcf
except:
    pindel = "/home/syngen/Software/pindel/pindel"
    pindel2vcf = "/home/syngen/Software/pindel/pindel2vcf"

import os
import sys
import datetime
import subprocess
import glob

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
        ori_info = Type + "\t" + Var_Type + "\t" + Interval + "\t" + Raw + "\t" + New + "\t" + fasta + "\t" + str(multiple) + "\t" + "Pindel"

        mes = ""
        for fn_line in fn:
            if "#" not in fn_line:
                pos = fn_line.split("\t")[1]
                ref = fn_line.split("\t")[3]
                alt = fn_line.split("\t")[7].split(";")[2]
                var_type = fn_line.split("\t")[7].split(";")[3]
                alt_reads = int(fn_line.split(":")[-1].split(",")[1])
                total_reads = int(fn_line.split(":")[-1].split(",")[1]) + int(fn_line.split(":")[-1].split(",")[0]) 
                ratio = float(alt_reads / total_reads)
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

        if not mes:
            mes = ori_info + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "\n"
    
    print(mes)

    return mes

def pindel_testing(fa, bam, info_record, *args):
    
    print("============= begin to call SNV using Pindel ==============")

    bam_config = ".".join(bam.split(".")[:-1]) + "_bam_config.txt"
    cmd = 'echo "%s\t300\tsample1\n" > %s' % (bam, bam_config)
    subprocess.call(cmd, shell=True)

    pindel_preout = ".".join(bam.split(".")[:-1]) + ".ipindel"
    cmd1 = "%s -f %s -i %s -c ALL -o %s -x 4 -w 1 -l false -M 10" % (pindel, fa, bam_config, pindel_preout)
    subprocess.call(cmd1, shell=True)

    specdate = datetime.datetime.now()
    date = str(specdate.year) + "-" + str(specdate.month) + "-" + str(specdate.day)
    vcf = ".".join(bam.split(".")[:-1]) + "_pindel.vcf"
    fa_prefix = fa.split("/")[2].split(".")[0]
    cmd2 = "%s -P %s -r %s -R %s -e 10 -co 10 -v %s -d %s" % (pindel2vcf, pindel_preout, fa, fa_prefix, vcf, date)
    print(cmd2)
    subprocess.call(cmd2, shell=True)
    mes = merge_files(info_record, vcf)

    ######### delete files we don't need #########
    os.remove(bam_config)
    for del_file in glob.glob("./data/*pindel_*"):
        os.remove(del_file)
    
    print("=================== Ending of Analysis  =====================")

    return mes

