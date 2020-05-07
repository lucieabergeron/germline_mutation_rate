# -*- coding: utf-8 -*-
"""
This script look for callability per trio:
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *
import pandas as pd

# Directories:
direct_denovo="{}/{}/de_novo_mutation/".format(path, sp)
direct_handling="{}/{}/vcf_handling/".format(path, sp)
direct="{}/{}/".format(path, sp)

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

# Import depth
depth_mean = pd.read_csv('{}depth.txt'.format(direct),sep=' ', index_col=None, header=None)

# The functions:
def call(direct_denovo, direct_handling, name, GQ_lim, DP_min, DP_max):
    """BCF tools for filters"""
    filt_cmd = "bcftools view -i 'MIN(FMT/DP)>={} && MIN(FMT/GQ)>={} && MAX(FMT/DP)<={}' {}inter_vcf/back_combine_genomicDBI_{}.g.vcf >> {}callability_filt_{}.g.vcf".format(DP_min, GQ_lim, DP_max, direct_handling, name, direct_denovo, name)
    """BCF tools for HomRef"""
    hom_cmd = "bcftools view -i 'FORMAT/PL[0:0]=0 && FORMAT/PL[1:0]=0' {}callability_filt_{}.g.vcf >> {}callability_filt_HomRef_{}.g.vcf".format(direct_denovo, name, direct_denovo, name)
    """Create a .sh files with the BCF functions and grep the number of sites."""
    file = open('{}callability_{}.sh'.format(direct_denovo, name),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --mem 64G \n')
    file.write('#SBATCH -c 16 \n')
    file.write('#SBATCH --time=20:00:00 \n')
    file.write('## BCFtools\n')
    file.write(filt_cmd)
    file.write('\n')
    file.write(hom_cmd)
    file.write('\n')
    file.write('## Grep: \n')
    file.write("echo Number with DP and GQ \n")
    file.write("grep -v '#' {}callability_filt_{}.g.vcf | wc -l \n".format(direct_denovo, name))
    file.write("echo Number of HomRef \n")
    file.write("grep -v '#' {}callability_filt_HomRef_{}.g.vcf | wc -l \n".format(direct_denovo, name))
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}callability_{}.out {}callability_{}.sh".format(direct_denovo, name, direct_denovo, name)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# Find the number of mutation per trios:
for name in trio_dir:
    mean_dp=depth_mean.loc[depth_mean[0] ==name][1]
    DP_min_calc=int(round(eval(DP_min)))
    DP_max_calc=int(round(eval(DP_max)))
    call(direct_denovo=direct_denovo, direct_handling=direct_handling, name=name, GQ_lim=GQ_lim, DP_min=DP_min_calc, DP_max=DP_max_calc)
    print("Callability for trio {}".format(name))
