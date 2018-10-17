#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import subprocess

tools_testing = "/home/syngen/zhangxu/tools_test/bin/tools_testing.py"

def exe_script(bps, ratios, dup_bps, dup_cps):

    for ratio in ratios:
        dir_name = os.getcwd() + "/data_" + str(ratio)
        os.mkdir(dir_name)
        os.chdir(dir_name)

        for bp in bps:
            res_file = "./data/total_info_" + str(bp) + "bp.txt"
            cmd = "%s -b %s -t 20 -d 10 --del %s --ins %s --inv %s --rep %s %s" % (tools_testing, bp, ratio, ratio, ratio, ratio, res_file)
            subprocess.call(cmd, shell=True)

        ###########################################
        for dup_bp in dup_bps:
            for dup_cp in dup_cps:
                dup_res_file = "./data/total_info_" + str(dup_bp) + "bp" + str(dup_cp) + "cp.txt"
                cmd2 = "%s -b %s -t 20 -p %s -d 10 --dup %s %s" % (tools_testing, dup_bp, dup_cp, ratio, dup_res_file)
                subprocess.call(cmd2, shell=True)
       ############################################

        new_dir = "/".join(dir_name.split("/")[:-1])
        os.chdir(new_dir)

if __name__ == "__main__":
    
    bps = [1, 10, 100, 200, 1000]
    ratios = [0.5, 0.1, 0.05, 0.01]

    dup_bps = [20, 60, 100]
    dup_cps = [2, 3]

    exe_script(bps, ratios, dup_bps, dup_cps)
