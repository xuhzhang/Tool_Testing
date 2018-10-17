#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import subprocess

tools_testing = "/home/syngen/zhangxu/tools_test/bin/tools_testing.py"

def exe_script(bps, ratios):

    for ratio in ratios:
        dir_name = os.getcwd() + "/data_" + str(ratio)
        os.mkdir(dir_name)
        os.chdir(dir_name)

        for bp in bps:
            res_file = "./data/total_info_" + str(bp) + "bp.txt"
            cmd = "%s -b %s -d 1 --del %s %s" % (tools_testing, bp, ratio, res_file)
            subprocess.call(cmd, shell=True)

        ###########################################
        #for dup_bp in dup_bps:
        #    dup_res_file = "./data/total_info_" + str(bp) + "bp" + str(cp) + "cp.txt"

         #   cmd2 = "%s -b %s -t 30 -p %s -d 10 --dup %s %s" % (tools_testing, dup_bp, dup_cp, ratio, dup_res_file)
         #   subprocess.call(cmd2, shell=True)
       ############################################

        new_dir = "/".join(dir_name.split("/")[:-1])
        os.chdir(new_dir)

if __name__ == "__main__":
    
    #bps = [1, 10, 100, 150, 200, 1000]
    #ratios = [0.5, 0.1, 0.05, 0.01]

    #dup_bps = [20, 30, 50]
    #dup_cps = [2, 3]

    bps = [3]
    ratios = [0.5]

    exe_script(bps, ratios)
