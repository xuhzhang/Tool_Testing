#!/bin/ust/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

try:
    from config import picard, samtools, gatk
except:
    picard = "/home/syngen/Software/picard.jar"
    samtools = "/home/syngen/Software/samtools-1.9/samtools"
    gatk = "/home/syngen/Software/gatk-4.0.8.1/gatk"

import os
import sys
import subprocess

def merge_files(info_record, hvcf, mvcf):
    
    ori_loc = info_record['info'][1].split('-')[0]
    Type = "Duplication"
    Var_Type = info_record['info'][0]
    Interval = info_record['info'][1]
    Raw = info_record['info'][2]
    New = info_record['info'][3]
    fasta = info_record['raw_fasta']
    multiple = info_record['multiple']
    ori_info = Type + "\t" + Var_Type + "\t" + Interval + "\t" + Raw + "\t" + New + "\t" + fasta + "\t" + str(multiple) + "\t" + "GATK"

    mes = ""
    with open(hvcf, "r") as fh:
        for hline in fh:
            if "#" not in hline:
                pos = hline.strip().split("\t")[1]
                ref = hline.strip().split("\t")[3]
                alt = hline.strip().split("\t")[4]
                var_type = "HaplotypeCaller"
                ratio = hline.strip().split("\t")[7].split(";")[1].split("=")[1]
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

    with open(mvcf, 'r') as fs:
        for sline in fs:
            if "#" not in sline:
                pos = sline.strip().split("\t")[1]
                ref = sline.strip().split("\t")[3]
                alt = sline.strip().split("\t")[4]
                var_type = "Mutect2"
                ratio = sline.strip().split("\t")[-1].split(":")[2]
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

    if not mes:
        mes = ori_info + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\n"

    return mes

def GATK_testing(fa, bam, info_record, *args):

    fa_dict = ".".join(fa.split(".")[:-1]) + ".dict"
    cmd1 = "java -jar %s CreateSequenceDictionary R=%s O=%s > /dev/null 2>&1" % (picard, fa, fa_dict)
    subprocess.call(cmd1, shell=True)
    cmd2 = "%s faidx %s > /dev/null 2>&1" % (samtools, fa)
    subprocess.call(cmd2, shell=True)

    hvcf = ".".join(bam.split(".")[:-1]) + "_HaplotypeCaller.vcf"
    cmd3 = "%s HaplotypeCaller -R %s -I %s -O %s > /dev/null 2>&1" % (gatk, fa, bam, hvcf)
    subprocess.call(cmd3, shell=True)

    cmd4 = "%s view -H %s 2>/dev/null | grep '^@RG' | awk -F ':' '{print $5}' > tmp.txt" % (samtools, bam)
    subprocess.call(cmd4, shell=True)
    
    with open('tmp.txt', 'r') as fr:
        bam_name = fr.read().strip()

    os.remove('tmp.txt')

    mvcf = ".".join(bam.split(".")[:-1]) + "_Mutect2.vcf"
    cmd5 = "%s Mutect2 -R %s -I %s -tumor %s -O %s > /dev/null 2>&1" % (gatk, fa, bam, bam_name, mvcf)
    subprocess.call(cmd5, shell=True)

    mes = merge_files(info_record, hvcf, mvcf)

    hvcf_index = hvcf + ".idx"
    mvcf_index = mvcf + ".idx"
    os.remove(hvcf_index)
    os.remove(mvcf_index)

    return mes

