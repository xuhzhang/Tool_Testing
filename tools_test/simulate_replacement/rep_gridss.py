#!/bin/usr/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

try:
    from config import gridss
except:
    gridss = "/home/syngen/Software/gridss/gridss-1.9.0-gridss-jar-with-dependencies.jar"

import os
import sys
import shutil
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
        ori_info = Type + "\t" + Var_Type + "\t" + Interval + "\t" + Raw + "\t" + New + "\t" + fasta + "\t" + str(multiple) + "\t" + "Gridss"

        mes = ""
        for fn_line in fn:
            if "#" not in fn_line and "LOW_QUAL" not in fn_line.split("\t")[6]:
                pos = fn_line.strip().split("\t")[1]
                ref = fn_line.strip().split("\t")[3]
                alt = fn_line.strip().split("\t")[4]
                var_type = fn_line.strip().split("\t")[7].split(";")[-2]
                ratio = "-"
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

        if not mes:
            mes = ori_info + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "\n"

    print(mes)

    return mes
                
def gridss_testing(fa, bam, info_record, *args):

    print("============= begin to call SNV using Delly ==============")

    vcf = ".".join(bam.split(".")[:-1]) + "_gridss.vcf"
    assemble_bam = ".".join(bam.split(".")[:-1]) + "_gridss.bam" 
    cmd = "java -jar %s REFERENCE_SEQUENCE=%s INPUT=%s OUTPUT=%s ASSEMBLY=%s" % (gridss, fa, bam, vcf, assemble_bam)
    subprocess.call(cmd, shell=True)

    mes = merge_files(info_record, vcf)

    print("=================== Ending of Analysis  =====================")

    bam_gridss = bam + ".gridss.working"
    assemble_gridss = assemble_bam + ".gridss.working"
    vcf_gridss = vcf + ".gridss.working"
    vcf_index = vcf + ".idx"
    del_file = [assemble_bam, vcf_index]
    del_dirs = [bam_gridss, vcf_gridss, assemble_gridss]

    for di in del_dirs:
        shutil.rmtree(di)

    for fi in del_file:
        os.remove(fi)

    return mes

