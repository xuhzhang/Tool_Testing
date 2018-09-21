#!/bin/usr/env python3
#! -*- coding:utf-8 -*-

__author__ = "zhangxu"

__author__ = "zhangxu"

try:
    from config import varscan, samtools
except:
    varscan = "/home/syngen/Software/varscan/VarScan.v2.3.9.jar"
    samtools = "/usr/local/bin/samtools"

import os
import sys
import subprocess

def merge_files(info_record, indel_vcf, snp_vcf):

    ori_loc = info_record['info'][1].split('-')[0]
    Type = "Inversion"
    Var_Type = info_record['info'][0]
    Interval = info_record['info'][1]
    Raw = info_record['info'][2]
    New = info_record['info'][3]
    fasta = info_record['raw_fasta']
    multiple = info_record['multiple']
    ori_info = Type + "\t" + Var_Type + "\t" + Interval + "\t" + Raw + "\t" + New + "\t" + fasta + "\t" + str(multiple) + "\t" + "VarScan"

    mes = ""
    with open(indel_vcf, "r") as f_ind:
        for iline in f_ind:
            if "Chrom" not in iline:
                pos = iline.strip().split("\t")[1]
                ref = iline.strip().split("\t")[2]
                alt = iline.strip().split("\t")[3]
                var_type = "indel"
                ratio = iline.strip().split("\t")[6]
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

    with open(snp_vcf, 'r') as f_snp:
        for sline in f_snp:
            if "Chrom" not in sline:
                pos = sline.strip().split("\t")[1]
                ref = sline.strip().split("\t")[2]
                alt = sline.strip().split("\t")[3]
                var_type = "snp"
                ratio = sline.strip().split("\t")[6]
                mes += ori_info + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + var_type + "\t" + str(ratio) + "\n"

    if not mes:
        mes = ori_info + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "-" + "\t" + "\n"

        print(mes)

    return mes

def varscan_testing(fa, bam, info_record, freq):

    print("============= begin to call SNV using VarScan2 ==============")

    varscan_indel_vcf = ".".join(bam.split(".")[:-1]) + "_varscan_indel.vcf"
    cmd = "samtools mpileup -f %s %s | java -jar %s pileup2indel --min-var-freq %f 1>%s" % (fa, bam, varscan, float(freq), varscan_indel_vcf)

    print(cmd)
    subprocess.call(cmd, shell=True)

    varscan_snp_vcf = ".".join(bam.split(".")[:-1]) + "_varscan_snp.vcf"
    cmd = "samtools mpileup -f %s %s | java -jar %s pileup2snp --min-var-freq %f 1>%s" % (fa, bam, varscan, float(freq), varscan_snp_vcf)

    print(cmd)
    subprocess.call(cmd, shell=True)

    mes = merge_files(info_record, varscan_indel_vcf, varscan_snp_vcf)

    print("=================== Ending of Analysis  =====================")

    return mes
