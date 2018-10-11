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
            cmd = "%s -b %s -t 30 -d 10 --del %s --ins %s --inv %s --rep %s" % (tools_testing, bp, ratio, ratio, ratio, ratio)
            subprocess.call(cmd, shell=True)

        new_dir = "/".join(dir_name.split("/")[:-1])
        os.chdir(new_dir)

if __name__ == "__main__":
    
    bps = [1, 10, 100, 150, 200, 1000]
    ratios = [0.01, 0.005, 0.001]

    exe_script(bps, ratios)
