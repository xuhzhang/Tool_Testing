#!/usr/bin/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

try:
    from config import bwa,samtools,samblaster,sambamba
except:
    bwa = "/usr/local/bin/bwa"
    samtools = "/usr/local/bin/samtools"
    samblaster = "/Users/xu/anaconda3/bin/samblaster"
    sambamba = "/Users/xu/anaconda3/bin/sambamba"

import subprocess

def mapping(ref, fq1, fq2):

    cmd1 = "%s index %s" % (bwa, ref)
    subprocess.call(cmd1, shell=True)

    cmd2 = "%s faidx %s" % (samtools,ref)
    subprocess.call(cmd2, shell=True)

    bam_out = "_".join(fq1.split("_")[:-1]) + "_bwa.bam"
    cmd3 = "%s mem -M -v 1 -R '@RG\\tID:IDwx\\tLB:LBwx\\tSM:SMwx\\tPL:ILLUMINA' %s %s %s | %s -M -q | %s view -S -f bam -l 0 /dev/stdin | %s sort -o %s /dev/stdin" % (bwa, ref, fq1, fq2, samblaster, sambamba, sambamba, bam_out)
    subprocess.call(cmd3, shell=True)

    flag_out = bam_out + ".flagstat"
    cmd4 = "%s flagstat %s > %s" % (samtools, bam_out, flag_out)
    subprocess.call(cmd4, shell=True)

    stats_out = bam_out + ".stats"
    cmd5 = "%s stats %s > %s" % (samtools, bam_out, stats_out)
    subprocess.call(cmd5, shell=True)

    return bam_out, flag_out, stats_out
