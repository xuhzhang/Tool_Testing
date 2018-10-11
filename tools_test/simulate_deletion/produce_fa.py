#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import random

def produce_fasta(prefix):
    
    if not os.path.exists("./data"):
        os.mkdir("./data")

    fasta_name = "./data/" + prefix + ".fa"

    if not os.path.exists(fasta_name):
        with open(fasta_name, 'w') as fw:
            fw.write(">ref\n")
            seq = ["A", "T", "G", "C"]
            for i in range(40000):
                random_base = random.choice(seq)
                fw.write(random_base)

    return fasta_name
