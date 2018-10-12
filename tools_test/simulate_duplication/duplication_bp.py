#!/usr/bin/env python3
# -*- coding:utf-8 -*-

__author__ = "zhangxu"

import random
import os
import json

def raw_seq(raw_fasta):

    with open(raw_fasta, 'r') as fr:
        fr.readline()
        info = fr.readline()

    return info

def duplicated_seq(raw_base, copy_number):

    copy_number = int(copy_number)
    new_base = ''
    for i in range(copy_number):
        new_base += raw_base

    return new_base

def duplication_bp(bp, raw_fasta, copy_number):

    bp = int(bp)
    info = raw_seq(raw_fasta)

    ###### write the duplicated bases to the new fasta ######
    new_fasta = "/".join(raw_fasta.split("/")[:-1]) + "/new_" + "_".join(raw_fasta.split("/")[-1].split("_")[1:])
    json_file = ".".join(new_fasta.split(".")[:-1]) + ".json"

    ### whether to produce new fasta or not ####
    if os.path.exists(json_file):
        with open(json_file, 'r') as fjr:
            info_record = json.load(fjr)
    else:
        ref = ">duplication_" + str(bp) + "bp.ref" + "\n"
        pos_bp = random.randint(0, (40000 - int(bp)))
        raw_base = info[pos_bp : (pos_bp+bp)]
        new_base = duplicated_seq(raw_base, copy_number)

        new_info = info[:pos_bp] + new_base + info[(pos_bp+bp):]
        with open(new_fasta, 'w') as fw:
            fw.write(ref)
            fw.write(new_info)
        ###### if the length of base is more than 10 bp ######
        ###### then it will be replaced by "-" ###############
        if int(bp) > 10:
            raw_base = "-"   
        if int(bp) * int(copy_number) > 30:
            new_base = "-"
        ###### use the dictionary to save the variated information ######
        info_record = {}
        var_interval = str(pos_bp+1) + "-" + str(pos_bp + bp)
        var_type = str(bp) + "bp" + str(copy_number) + "cp"
        info_record['info'] = [var_type, var_interval, raw_base, new_base]
        info_record['fasta'] = new_fasta
    
        ###### add the record file to avoid info lost ########
        with open(json_file, 'w') as fjw:
            json.dump(info_record, fjw)

    return info_record
