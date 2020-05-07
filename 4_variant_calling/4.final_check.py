# -*- coding: utf-8 -*-
"""
This script checked if:
 the gather genotype and gather back_combined exist
 move to log_files the useless files
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
direct= "{}/{}/vcf_files/".format(path, sp)

# Dictionary made of tuples of sample name and vcf file:
f = open('{}/{}/vcf_files.txt'.format(path, sp))
vcf_dir = {}
for line in f:
    name = line.split('_')[0]
    vcf = line.split()[0]
    if name not in vcf_dir:
        vcf_dir[name] = []
    vcf_dir[name] = vcf


##################################################
# What you run  ##################################
##################################################

if os.path.exists("{}back_combine_genomicDBI_gather.g.vcf".format(direct)) and os.path.exists("{}genotype_genomicDBI_gather.g.vcf".format(direct)):
    print("Gathered back combine and genotype exist --> ready for vcf_handeling")
## Something else...
    mv_vcf = "mv {}*_chr*.g.vcf* {}inter_vcf/".format(direct, direct)
    subprocess.call(mv_vcf, shell=True)
    mv_vcf1 = "mv {}*_scaff*.g.vcf* {}inter_vcf/".format(direct, direct)
    subprocess.call(mv_vcf1, shell=True)
    print("Move the inter vcf files")
    mv_log = "mv {}gather_genomicDBImport* {}gather.log".format(direct, direct)
    subprocess.call(mv_log, shell=True)
    print("Move the gather log files")
else:
    if os.path.exists("{}genotype_genomicDBI_gather.g.vcf".format(direct)):
        print("Gathered genotype exists --> PROBLEM with back combine")
    elif os.path.exists("{}back_combine_genomicDBI_gather.g.vcf".format(direct)):
        print("Gathered back combine exists --> PROBLEM with genotype")

