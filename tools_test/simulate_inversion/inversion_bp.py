#!/usr/bin/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

import random

def raw_seq(raw_fasta):

    with open(raw_fasta, 'r') as fr:
        fr.readline()
        info = fr.readline()

    return info

def inversed_seq(raw_base):

    complement_seq = ''
    for seq in raw_base:
        if seq == 'A':
            complement_seq += 'T'
        elif seq == 'T':
            complement_seq += 'A'
        elif seq == 'C':
            complement_seq += 'G'
        elif seq == 'G':
            complement_seq += 'C'

    reversed_seq = ''
    for reseq in reversed(complement_seq):
        reversed_seq += reseq

    return reversed_seq

def inversion_bp(bp, raw_fasta):

    bp = int(bp)
    info = raw_seq(raw_fasta)

    ###### write the replaced bases to the new fasta ######
    new_fasta = "./data/inversion_" + str(bp) + "bp.fasta"
    ref = ">inversion_" + str(bp) + "bp.ref" + "\n"
    pos_bp = random.randint(0, (40000 - int(bp)))

    raw_base = info[pos_bp:(pos_bp+bp)]
    new_base = inversed_seq(raw_base)

    ##############################################################

    while raw_base == new_base:
        pos_bp = random.randint(0, 40000 - int(bp))
        raw_base = info[pos_bp:(pos_bp+bp)]
        new_base = inversed_seq(raw_base)

    ###### if the variated reads length is more than 10bp, #######
    ###### then raw_base and new_base will be presented as "-" ###
    #if int(bp) > 10:
    #    raw_base = "-"
    #    new_base = "-"
    #else:
    #    raw_base = info[pos_bp:(pos_bp+bp)]
    #    new_base = inversed_seq(raw_base)

    ##############################################################

    new_info = info[:pos_bp] + new_base + info[(pos_bp+bp):]
    with open(new_fasta, 'w') as fw:
        fw.write(ref)
        fw.write(new_info)

    ###### use the dictionary to save the variated information ######
    info_record = {}
    var_interval = str(pos_bp + 1) + "-" + str(pos_bp + bp)
    var_type = str(bp) + "bp"
    info_record['info'] = [var_type, var_interval, raw_base, new_base]
    info_record['fasta'] = new_fasta

    return info_record
