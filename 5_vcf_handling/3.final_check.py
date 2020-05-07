# -*- coding: utf-8 -*-
"""
This script do a final check of the file, move the log files and summarized the depth.

"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
direct = "{}/{}/vcf_handling/".format(path, sp)
dp_dir = "{}/{}/".format(path, sp)

# Dictionary:
f = open('{}/{}/pedigree.ped'.format(path, sp))
trio_dir = {}
for line in f:
    off = line.split()[1]
    fa = line.split()[2]
    mo = line.split()[3]
    name = off
    if name not in trio_dir:
        trio_dir[name] = []
    trio_dir[name].append((off, fa, mo))


##################################################
# What you run  ##################################
##################################################

# Check if all the last file exists and move the log:
geno_exist=[]
back_exist=[]
for name in trio_dir:
    geno_exist.append(os.path.exists("{}trio_{}_vcf_table.out".format(direct, name)))
    back_exist.append(os.path.exists("{}depth_{}.out".format(direct, name)))

if all(geno_exist) and all(back_exist):
    print("All the table and depth exist --> Ready for de novo detection")
    print("Move the log files")
    sub_cmd = "mv -t {}log_file {}*geno_handling* {}relatedness.* {}pedigree* {}DP* {}back_handling* {}*coverage*".format(direct, direct, direct, direct, direct, direct, direct)
    subprocess.call(sub_cmd, shell=True)
    mv_vcf = "mv -t {}inter_vcf/ {}*genotype_genomicDBI* {}*back_combine_genomicDBI*".format(direct, direct, direct)
    subprocess.call(mv_vcf, shell=True)
    print("Move the inter vcf files")
    print("Get the summary of each depth in {}depth.txt".format(dp_dir))
    file = open('{}depth.txt'.format(dp_dir),'w')
    file.close()
    for name in trio_dir:
        sum_cmd="echo \"{} $(grep 'tot' {}depth_{}.out | cut -d' ' -f3)\">>{}depth.txt".format(name, direct, name, dp_dir)
        subprocess.call(sum_cmd, shell=True)
else:
    print("PROBLEM --> Some of the tables file doesn't exist")
