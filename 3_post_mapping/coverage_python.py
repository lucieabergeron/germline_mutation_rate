### -*- coding: utf-8 -*-
##"""
##This script get the coverage.
##"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import pandas as pd

# Dictionary made of tuples of sample name and merged bamfile:
def average_column(csv, name):
    f = open(csv)
    # Basic stats:
    average = 0
    Sum = 0
    Max = 0
    Min = 10000
    row_count = 0
    # Sup a certain cov:
    sup_10 = 0
    sup_30 = 0
    sup_50 = 0
    sup_70 = 0
    for line in f:
        cov = line.split()[2]
        n=float(cov)
        Sum += n
        row_count += 1
        Max = max(Max,n)
        Min = min(Min,n)
        if n >= 70 :
            sup_70 += 1
            sup_50 += 1
            sup_30 += 1
            sup_10 += 1
        elif n >= 50:
            sup_50 += 1
            sup_30 += 1
            sup_10 += 1
        elif n >= 30:
            sup_30 += 1
            sup_10 += 1
        elif n >= 10:
            sup_10 += 1
    average = Sum / row_count
    average = round(average,2)
    summary = "{} {} {} {} . {} . {} . {}".format(average, Max, Min, sup_10, sup_30, sup_50, sup_70)
    cmd = "echo \"{} {} \" >> {}/{}/bam_files/coverage/summary_coverage.txt".format(name, summary, path, sp)
    subprocess.call(cmd, shell=True)
    f.close()

# Function call
