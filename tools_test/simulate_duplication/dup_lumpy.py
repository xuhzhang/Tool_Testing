#!/bin/usr/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

try:
    from config import lumpyexpress, samtools,extractSplitReads_BwaMem
except:
    lumpyexpress = "/home/syngen/Software/lumpy-sv/bin/lumpyexpress"
    samtools = "/home/syngen/Software/samtools-1.9/samtools"
    extractSplitReads_BwaMem = "/home/syngen/Software/lumpy-sv/scripts/extractSplitReads_BwaMem"

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
        ori_info = Type + "\t" + Var_Type + "\t" + Interval + "\t" + Raw + "\t" + New + "\t" + fasta + "\t" + str(multiple) + "\t" + "lumpy"

        mes = ""
        for fn_line in fn:
            if "#" not in fn_line:
                pos = fn_line.strip().split("\t")[1]
                ref = fn_line.strip().split("\t")[3]
                alt = fn_line.strip().split("\t")[7].split(";")[2]
                var_type = fn_line.strip().split("\t")[7].split(";")[0] 
                ratio = "-" 
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

        if not mes:
            mes = ori_info + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "\n"

    print(mes)

    return mes

def lumpy_testing(_fa, bam, info_record, *args):

    print("============= begin to call SV using lumpyexpress ==============")

    split_bam = ".".join(bam.split(".")[:-1]) + "_splitters.bam"
    cmd1 = "%s view -h %s | %s -i stdin | %s view -Sb - | %s sort > %s" % (samtools, bam, extractSplitReads_BwaMem, samtools, samtools, split_bam)
    subprocess.call(cmd1, shell=True)

    print(cmd1)

    disordants_bam = ".".join(bam.split(".")[:-1]) + "_disordants.bam"
    cmd2 = "%s view -b -F 1294 %s | %s sort > %s" % (samtools, bam, samtools, disordants_bam)
    subprocess.call(cmd2, shell=True)

    print(cmd2)

    vcf = ".".join(bam.split(".")[:-1]) + "_lumpy.vcf"
    cmd3 = "%s -B %s -S %s -D %s -o %s" % (lumpyexpress, bam, split_bam, disordants_bam, vcf)
    subprocess.call(cmd3, shell=True)

    print(cmd3)

    mes = merge_files(info_record, vcf)

    ### delete the unnecessary files ###
    os.remove(split_bam)
    os.remove(disordants_bam)

    print("=================== Ending of Analysis  =====================")

    return mes

