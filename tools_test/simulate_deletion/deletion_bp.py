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

def deletion_bp(bp, raw_fasta):

    info = raw_seq(raw_fasta)

    ###### write the deleted bases to the new fasta ######
    new_fasta = "/".join(raw_fasta.split("/")[:-1]) + "/new_" + "_".join(raw_fasta.split("/")[-1].split("_")[1:])
    json_file = ".".join(new_fasta.split(".")[:-1]) + ".json"

    ### whether to produce new fasta or not ####
    if os.path.exists(json_file):
        with open(json_file, 'r') as fjr:
            info_record = json.load(fjr)
    else:
        ref = ">del_" + str(bp) + "bp.ref" + "\n"
        pos_bp = random.randint(0, (40000 - int(bp)))
        raw_base = info[pos_bp : (pos_bp + int(bp))]
        new_info = info[:pos_bp] + info[(pos_bp + int(bp)):]
        with open(new_fasta, 'w') as fw:
            fw.write(ref)
            fw.write(new_info)
        ###### if the length of base is more than 10 bp ######
        ###### then it will be replaced by "-" ###############
        if int(bp) > 10:
            raw_base = "-"
        ###### use the dictionary to save the variated information ######
        info_record = {}
        var_interval = str(pos_bp+1) + "-" + str(pos_bp + int(bp))
        var_type = str(bp) + "bp"
        info_record["info"] = [var_type, var_interval, raw_base]
        info_record["fasta"] = new_fasta

        ###### add the record file to avoid info lost ########
        with open(json_file, 'w') as fjw:
            json.dump(info_record, fjw)

    return info_record
